# Annotate and visualize the white spruce mitochondrial genome
# Copyright 2014 Shaun Jackman

name=pg29mt-concat

# Cycas taitungensis mitochondrial DNA, complete genome
ref=NC_010303

all: $(name).gbk.png

clean:
	rm -f $(name).gb $(name).gbk $(name).gff $(name).png

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Fetch data from NCBI

resources/%.gbk:
	curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&val=$*" >$@

# asn,faa,ffn,fna,frn,gbk,gff,ptt,rnt,rpt,val
resources/%:
	mkdir -p resources
	curl -fsS ftp://ftp.ncbi.nih.gov/genomes/MITOCHONDRIA/Metazoa/$* >$@

%.faa: %.gbk
	bin/gbk-to-faa <$< >$@

%.frn: resources/%.frn
	sed 's/^>.*\[gene=/>/;s/\].*$$//' $< >$@

# RepeatModeler

%.nin: %.fa
	/usr/local/opt/repeatmodeler/BuildDatabase -name $* -engine ncbi $<

%.rm: %.nin
	RepeatModeler -database $*
	touch $@

RepeatModeler.fa: ThuApr31728402014/consensi.fa.classified
	cp -a $< $@

# MAKER

maker_bopts.ctl:
	maker -BOPTS

maker_exe.ctl:
	maker -EXE

rmlib.fa: PICEAGLAUCA_rpt2.0.fa RepeatModeler.fa
	cat $^ >$@

%.maker.output/stamp: maker_opts.ctl %.fa $(ref).frn $(ref).faa rmlib.fa
	maker -fix_nucleotides -cpus 4
	touch $@

%.gff: %.maker.output/stamp
	gff3_merge -s -g -n -d $*.maker.output/$*_master_datastore_index.log \
		|sed 's/trnascan-[^-]*-processed-gene/trn/g;s/-tRNA-1//g' \
		|sed '/rrn/s/mRNA/rRNA/;/trn/s/mRNA/tRNA/' >$@

%.gb: %.gff %.fa
	bin/gff_to_genbank.py $^ >$@

%.gbk: %-header.gbk %.gb
	(cat $< && sed -n '/^FEATURES/,$${s/Name=/gene=/;s/-gene//;p;}' $*.gb) >$@

# OrganellarGenomeDRAW

%.gbk.png: %.gbk
	drawgenemap --format png --infile $< --outfile $<
