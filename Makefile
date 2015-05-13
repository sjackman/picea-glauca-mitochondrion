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
	brew install aragorn bedtools edirect genometools maker prokka repeatmodeler trnascan
	pip install biopython seqmagick

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

# Copy local data

#PICEAGLAUCA_rpt2.0.fa: /genesis/extscratch/seqdev/PG/data/PICEAGLAUCA_rpt2.0
#	cp -a $< $@

# Fetch data from NCBI

cds_aa.orig.fa cds_na.orig.fa: %.fa:
	esearch -db nuccore -query $(edirect_query) \
		|efetch -format fasta_$* >$@

cds_aa.fa cds_na.fa: %.fa: %.orig.fa
	sed -E 's/^>(.*gene=([^]]*).*)$$/>\2|\1/' $< \
		|seqmagick -q convert \
		--pattern-exclude '^lcl|^orf|^ORF|hypothetical|putative|unnamed' \
		--pattern-replace '^apt' 'atp' \
		--pattern-replace '^coxIII' 'cox3' \
		--pattern-replace '^coxII' 'cox2' \
		--pattern-replace '^coxI' 'cox1' \
		--pattern-replace '^cytb' 'cob' \
		--pattern-replace '^nd' 'nad' \
		--pattern-replace '^yejU' 'ccmC' \
		--pattern-replace '^yejV' 'ccmB' \
		--deduplicate-taxa \
		- $@

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

# Barrnap

%.barrnap.gff: %.fa
	barrnap --kingdom mitochondria --threads $t $< >$@

# Annotate rRNA using RNAmmer

%.rnammer.gff2: %.fa
	mkdir -p rnammer
	rnammer -S bac -gff $@ -xml rnammer/$*.xml -f rnammer/$*.fa -h rnammer/$*.hmm $<

# Convert GFF2 to GFF3
%.rnammer.gff: %.rnammer.gff2
	bin/convert_RNAmmer_to_gff3.pl --input=$< >$@

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

# Remove mRNA records
%.nomrna.gff: %.gff
	gsed '/\tmRNA\t/d' $< >$@

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
		/\tCDS\t/{/gene=/!s/ID=[^_]*_([0-9]*)/Name=orf\1;&/;}; \
		' $< >$@

# Report the genes annotated by Prokka
prokka/%.gff.gene: prokka/%.gff
	ruby -we 'ARGF.each { |s| \
		puts $$1 if s =~ /\tgene\t.*gene=([^;]*)/ \
	}' $< >$@

# Convert to GenBank format
%.gb: %.nomrna.gff %.fa
	bin/gff_to_genbank.py $^
	mv $*.nomrna.gb $@

%.gbk: %-header.gbk %.gb
	(cat $< && sed -En '/^FEATURES/,$$ { \
		s/Name="(trn[^-]*).*"/gene="\1"/; \
		s/Name="([^|]*).*"/gene="\1"/; \
		p; }' $*.gb) >$@

# Merge MAKER and Prokka annotations using bedtools

%.gff: %.prokka.gff %.maker.gff
	bedtools intersect -v -header -a $< -b $*.maker.gff \
		|sed 's/CDS/mRNA/' \
		|gt gff3 -sort $*.maker.gff - >$@

# OrganellarGenomeDRAW

%.gbk.png: %.gbk
	drawgenemap --format png --infile $< --outfile $<

# Report the annotated genes

%.gff.gene: %.gff
	bin/gff-gene-name $< >$@

# Split the GenBank file into one sequence per file
gbk/%.00.gbk: %.gbk
	gcsplit -sz -f $*. --suppress-matched $< '/\/\//' '{*}'
	mkdir -p $(@D)
	for i in $*.[0-9][0-9]; do mv $$i gbk/$$i.gbk; done

# Combine the OGDraw images into a single image
%.gbk.montage.png: gbk/%.*.gbk.png
	montage -geometry 1200x600 $^ $@

# GenomeTools sketch

%.gff.png: %.gff
	gt sketch $@ $<

# Render HTML from RMarkdown
%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$<", output_format = "html_document")'

# Dependencies

genes.html: pg29mt-scaffolds.gff

repeats.html: pg29mt-scaffolds.maker.repeat.gff
