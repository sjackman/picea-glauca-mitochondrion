# Annotate and visualize the white spruce mitochondrial genome
# Copyright 2015 Shaun Jackman

# Name of the assembly
name=pg29mt-scaffolds

# Number of threads
t=4

# Green plant mitochondria
edirect_query='Viridiplantae[Organism] mitochondrion[Title] (complete genome[Title] OR complete sequence[Title])'

all: $(name).gff $(name).evidence.gff $(name).repeat.gff $(name).gff.gene $(name).gbk.png

clean:
	rm -f $(name).gb $(name).gbk $(name).gff $(name).png

install-deps:
	pip install seqmagick

.PHONY: all clean install-deps
.DELETE_ON_ERROR:
.SECONDARY:

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
	/usr/local/opt/repeatmodeler/BuildDatabase -name $* -engine ncbi $<

%.RepeatModeler.fa: %.nin
	RepeatModeler -database $*
	cp -a RM_*/consensi.fa.classified $@

# Barrnap

%.rrna.gff: %.fa
	barrnap --kingdom mitochondria --threads $t $< >$@

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

%.repeat.gff: %.maker.output/stamp
	cat `find $*.maker.output -name query.masked.gff` >$@

%.evidence.orig.gff: %.maker.output/stamp
	gff3_merge -s -n -d $*.maker.output/$*_master_datastore_index.log >$@

%.orig.gff: %.maker.output/stamp
	gff3_merge -s -g -n -d $*.maker.output/$*_master_datastore_index.log >$@

%.gff: %.orig.gff
	gsed -E 's/Name=trnascan-[^-]*-noncoding-([^-]*)-gene/Name=trn\1/g; \
		s/Name=([^;]*)S_rRNA/Name=rrn\1/g' \
		$^ >$@

# Add the rRNA annotations to the GFF file
$(name).gff: $(name).rrna.gff

# Remove mRNA records
%.nomrna.gff: %.gff
	gsed '/\tmRNA\t/d' $< >$@

# Convert to GenBank format
%.gb: %.nomrna.gff %.fa
	bin/gff_to_genbank.py $^
	mv $*.nomrna.gb $@

%.gbk: %-header.gbk %.gb
	(cat $< && sed -En '/^FEATURES/,$$ { \
		s/Name="(trn[^-]*).*"/gene="\1"/; \
		s/Name="([^|]*).*"/gene="\1"/; \
		p; }' $*.gb) >$@

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

# Render HTML from RMarkdown
%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$<", output_format = "html_document")'

# Dependencies

genes.html: pg29mt-scaffolds.gff

repeats.html: pg29mt-scaffolds.repeat.gff
