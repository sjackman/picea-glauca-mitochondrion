# Annotate and visualize the white spruce mitochondrial genome
# Copyright 2015 Shaun Jackman

# Name of the assembly
name=pg29mt-scaffolds

# Number of threads
t=4

# Green plant mitochondria
edirect_query='Viridiplantae[Organism] mitochondrion[Title] (complete genome[Title] OR complete sequence[Title])'

all: $(name).gff $(name).gbk $(name).gbk.png \
	$(name).maker.evidence.gff $(name).maker.repeat.gff \
	$(name).maker.gff.gene $(name).prokka.gff.gene $(name).gff.gene \
	genes.html repeats.html

clean:
	rm -f $(name).gb $(name).gbk $(name).gff $(name).png

install-deps:
	brew install aragorn bedtools edirect genometools gnu-sed maker ogdraw prokka repeatmodeler rnammer trnascan
	pip install bcbio-gff biopython seqmagick

.PHONY: all clean install-deps
.DELETE_ON_ERROR:
.SECONDARY:

# Download scripts

bin/convert_RNAmmer_to_gff3.pl:
	curl -o $@ https://raw.githubusercontent.com/jorvis/biocode/master/gff/convert_RNAmmer_to_gff3.pl
	chmod +x $@

bin/aragorn_out_to_gff3.py:
	curl -o $@ https://raw.githubusercontent.com/bgruening/galaxytools/master/tools/rna_tools/trna_prediction/aragorn_out_to_gff3.py
	chmod +x $@

# BLAST

# Align the scaffolds to the nt database
%.blastn: %.fa
	blastn -db nt -query $< -out $@

# Copy local data

#PICEAGLAUCA_rpt2.0.fa: /genesis/extscratch/seqdev/PG/data/PICEAGLAUCA_rpt2.0
#	cp -a $< $@

# Fetch data from NCBI

cds_aa.orig.fa cds_na.orig.fa: %.fa:
	esearch -db nuccore -query $(edirect_query) \
		|efetch -format fasta_$* >$@

cds_aa.fa cds_na.fa: %.fa: %.orig.fa
	seqmagick -q convert \
		--pattern-exclude 'gene=orf|hypothetical|putative|unnamed' \
		--pattern-replace 'gene=apt' 'gene=atp' \
		--pattern-replace 'gene=ccmFn' 'gene=ccmFN' \
		--pattern-replace 'gene=coxIII' 'gene=cox3' \
		--pattern-replace 'gene=coxII' 'gene=cox2' \
		--pattern-replace 'gene=coxI' 'gene=cox1' \
		--pattern-replace 'gene=cytb' 'gene=cob' \
		--pattern-replace 'gene=nd' 'gene=nad' \
		--pattern-replace 'gene=yejU' 'gene=ccmC' \
		--pattern-replace 'gene=yejV' 'gene=ccmB' \
		--pattern-replace 'gene=18S rRNA' 'gene=40' \
		$< - \
	|gsed -E '/protein=[^]]*intron[^]]*ORF/s/gene=/gene=ymf/; \
		s/^>(.*gene=([^]]*).*)$$/>\2|\1/' \
	|seqmagick -q convert --pattern-exclude '^lcl' --deduplicate-taxa - $@

# Cycas taitungensis
NC_010303.1.json: %.json:
	bionode-ncbi search nuccore $* >$@

%.uid: %.json
	json uid <$< >$@

%.fa: %.uid
	curl http://togows.org/entry/nucleotide/`<$<`.fasta |seqtk seq >$@

%.gb: %.uid
	curl -o $@ http://togows.org/entry/nucleotide/`<$<`.gb

%.gff: %.uid
	curl -o $@ http://togows.org/entry/nucleotide/`<$<`.gff

# Organelle Genome Resources
# http://www.ncbi.nlm.nih.gov/genome/organelle/
mitochondrion/all: mitochondrion/mitochondrion.1.1.genomic.fna.gz mitochondrion/mitochondrion.1.genomic.gbff.gz mitochondrion/mitochondrion.1.protein.faa.gz mitochondrion/mitochondrion.1.protein.gpff.gz mitochondrion/mitochondrion.1.rna.fna.gz mitochondrion/mitochondrion.1.rna.gbff.gz

mitochondrion/mitochondrion.%.gz:
	mkdir -p $(@D)
	curl -o $@ ftp://ftp.ncbi.nlm.nih.gov/refseq/release/$@

# RepeatModeler

%.nin: %.fa
	BuildDatabase -name $* -engine ncbi $<

%.RepeatModeler.fa: %.nin
	RepeatModeler -database $*
	cp -a RM_*/consensi.fa.classified $@

# ARAGORN

# Annotate tRNA using ARAGORN and output TSV
%.aragorn.tsv: %.fa
	aragorn -gcstd -i -l -w -o $@ $<

# Annotate tRNA using ARAGORN and output text
%.aragorn.txt: %.fa
	aragorn -gcstd -i -l -o $@ $<

# Convert ARAGORN output to GFF
%.aragorn.gff: %.aragorn.tsv
	bin/aragorn_out_to_gff3.py --full <$< |gt gff3 -sort >$@

# Annotate tRNA using tRNAscan-SE
%.trnascan.orig.tsv: %.fa
	tRNAscan-SE -O -o $@ -f $*.trnascan.txt $<

# Barrnap

%.barrnap.gff: %.fa
	barrnap --kingdom mitochondria --threads $t $< >$@

# Annotate rRNA using RNAmmer

%.rnammer.gff2: %.fa
	mkdir -p rnammer
	rnammer -S bac -gff $@ -xml rnammer/$*.xml -f rnammer/$*.fa -h rnammer/$*.hmm $<

# Convert GFF2 to GFF3
%.rnammer.gff: %.rnammer.gff2
	bin/convert_RNAmmer_to_gff3.pl --input=$< \
		|sed -E 's/ID=([^s]*)s_rRNA_([0-9]*)/Name=rrn\1;&/' >$@

# MAKER

maker_bopts.ctl:
	maker -BOPTS

maker_exe.ctl:
	maker -EXE

rmlib.fa: PICEAGLAUCA_rpt2.0.fa $(name).RepeatModeler.fa
	cat $^ >$@

%.maker.output/stamp: maker_opts.ctl %.fa cds_aa.fa rmlib.fa
	maker -fix_nucleotides -cpus $t
	touch $@

%.maker.repeat.gff: %.maker.output/stamp
	cat `find $*.maker.output -name query.masked.gff` >$@

%.maker.evidence.gff: %.maker.output/stamp
	gff3_merge -s -n -d $*.maker.output/$*_master_datastore_index.log >$@

%.maker.orig.gff: %.maker.output/stamp
	gff3_merge -s -g -n -d $*.maker.output/$*_master_datastore_index.log >$@

%.maker.gff: %.maker.orig.gff
	gt gff3 -sort $^ \
	|gsed -E ' \
		/\tintron\t/d; \
		s/Name=trnascan-[^-]*-noncoding-([^-]*)-gene/Name=trn\1/g; \
		/\trRNA\t/s/ID=([^;]*)s_rRNA/Name=rrn\1;&/g' \
	|gt gff3 -addintrons -sort - >$@

# Add the rRNA annotations to the GFF file
$(name).maker.gff: $(name).rnammer.gff

# Add the tRNA annotations to the GFF file
$(name).maker.gff: $(name).aragorn.gff

# Prokka

# Convert the FASTA file to the Prokka FASTA format
cds_aa.prokka.fa: %.prokka.fa: %.fa
	sed -E 's/^>([^ ]*) .*gene=([^]]*).*protein=([^]]*).*$$/>\1 ~~~\2~~~\3/; \
		s/^-//' $< >$@

# Annotate genes using Prokka
prokka/%.gff: %.fa cds_aa.prokka.fa
	prokka --kingdom bac --gcode 1 --addgenes --proteins cds_aa.prokka.fa --rnammer \
		--cpus $t \
		--genus Picea --species 'glauca mitochondrion' \
		--locustag OU3MT \
		--force --outdir prokka --prefix $* \
		$<

# Remove the FASTA section from the Prokka GFF file
%.prokka.gff: prokka/%.gff
	gsed -E '/^##FASTA/,$$d; \
		s/gene=([^;]*)/Name=\1;&/; \
		/\tgene\t/{/gene=/!s/ID=[^_]*_([0-9]*)/Name=orf\1;&/;}; \
		/\tCDS\t/{/gene=/!s/ID=[^_]*_([0-9]*)/Name=orf\1;&/; \
			s/CDS/mRNA/;p; \
			s/mRNA/CDS/;s/Parent=[^;]*;//;s/ID=/Parent=/;}; \
		' $< >$@

# Report the genes annotated by Prokka
prokka/%.gff.gene: prokka/%.gff
	ruby -we 'ARGF.each { |s| \
		puts $$1 if s =~ /\tgene\t.*gene=([^;]*)/ \
	}' $< >$@

# Remove mRNA and ORF annotations before converting to GBK
%.pregbk.gff: %.gff
	gsed -E '/\tmRNA\t|Name=orf/d;/^[0-9]\t/s/^/0/' $< |uniq >$@

# Add leading zeros to the FASTA IDs
%.pregbk.fa: %.fa
	sed '/^>[0-9]$$/s/^>/>0/' $< >$@

# Convert to GenBank format
%.gb: %.pregbk.gff %.pregbk.fa
	bin/gff_to_genbank.py $^
	sed -e '/DEFINITION/{h;s/$$/ mitochondrion/;}' \
		-e '/ORGANISM/{g;s/DEFINITION/  ORGANISM/;}' \
		$*.pregbk.gb >$@

%.gbk: %-header.gbk %.gb
	(cat $< && sed -En '/^FEATURES/,$$ { \
		s/Name="([^|]*).*"/gene="\1"/; \
		p; }' $*.gb) >$@

# Merge MAKER and Prokka annotations using bedtools

%.gff: %.prokka.gff %.maker.gff
	bedtools intersect -v -header -a $< -b $*.maker.gff \
		|gt gff3 -sort $*.maker.gff - >$@

# OrganellarGenomeDRAW

%.gbk.png: %.gbk
	drawgenemap --density 150 --format png --infile $< --outfile $<

# Report the annotated genes

# Extract the names of genes from a GFF file
%.gff.gene: %.gff
	bin/gff-gene-name $< >$@

# Extract DNA sequences of GFF gene features from a FASTA file
%.gff.gene.fa: %.gff %.fa
	gt extractfeat -type gene -matchdescstart -retainids -seqid -seqfile $*.fa $< >$@

# Extract DNA sequences of GFF CDS features from a FASTA file
%.gff.CDS.fa: %.gff %.fa
	gt extractfeat -type CDS -join -matchdescstart -retainids -seqid -seqfile $*.fa $< >$@

# Translate protein sequences of GFF CDS features from a FASTA file
%.aa.fa: %.CDS.fa
	gt seqtranslate -reverse no -fastawidth 0 $< \
	|sed -n '/ (1+)$$/{s/ (1+)$$//;p;n;p;n;n;n;n;}' >$@

# UniqTag

# Generate UniqTag from DNA or amino acid sequence
%.uniqtag: %.fa
	uniqtag $< >$@

# Extract sequences of GFF intron features
%.gff.intron.fa: %.gff %.fa
	gt extractfeat -type intron -matchdescstart -retainids -seqid -seqfile $*.fa $< >$@

# Extract sequences of GFF intron features plus flanking sequence
%.gff.intron.flank100.fa: %.gff %.fa.fai
	(awk '$$3 == "intron"' $<; \
		awk '$$3 == "intron"' $< |bedtools flank -b 100 -i stdin -g $*.fa.fai) \
	|sort -k1,1n -k4,4n - \
	|bedtools merge -i stdin \
	|bedtools getfasta -bed stdin -fi $*.fa -fo $@

# Split the GenBank file into one sequence per file
gbk/%.00.gbk: %.gbk
	gcsplit -sz -f $*. --suppress-matched $< '/\/\//' '{*}'
	mkdir -p $(@D)
	for i in $*.[0-9][0-9]; do mv $$i gbk/$$i.gbk; done

# Combine the OGDraw images into a single image
%.gbk.montage.png: gbk/%.*.gbk.png
	montage -tile 3 -geometry +0+0 -units PixelsPerInch -density 1200 $^ gbk/$*.00.gbk_legend.png $@

# GenomeTools sketch
%.gff.png: %.gff
	gt sketch $@ $<

# Convert GFF to GTF
%.gtf: %.gff
	gt -q gff3_to_gtf $< >$@

# Render HTML from RMarkdown
%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<", output_format = "html_document")'
	mogrify -units PixelsPerInch -density 300 $*_files/figure-html/*.png

# Dependencies

genes.html: pg29mt-scaffolds.gff

repeats.html: pg29mt-scaffolds.maker.repeat.gff
