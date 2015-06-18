# Written by Shaun Jackman @sjackman

# For shell redirection <(gunzip -c reads.fq.gz)
SHELL=/bin/bash

all: abundance.html

install-deps:
	brew install kallisto samtools

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:

# Create symlinks

pg29mt-scaffolds.gff:
	ln -sf ../$@ .

# Convert BAM to FASTQ

%_1.fq.gz: %.bam
	time samtools view -uf64 $< |samtools bam2fq - |gzip >$@

%_2.fq.gz: %.bam
	time samtools view -uf128 $< |samtools bam2fq - |gzip >$@

# Kallisto

# Index the reference
%.idx: %.fa
	kallisto index -i $@ $<

# Quantify the transcripts
%/abundance.txt: transcripts.idx %_1.fq.gz %_2.fq.gz
	time kallisto quant -o $(@D) -i $^

# Salmon

# Index the reference
%/bwaidx.bwt: %.fa
	salmon index -t $< -i $(@D)

# Quantify the transcripts
%_quant/quant.sf: transcripts/bwaidx.bwt %_1.fq.gz %_2.fq.gz
	salmon quant -i $(<D) -l IU -1 <(gunzip -c $*_1.fq.gz) -2 <(gunzip -c $*_2.fq.gz) -o $(@D)

# RMarkdown

# Render the report to HTML
%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<", output_format = "html_document")'

# Dependencies

abundance.html: \
	tissues.tsv \
	Bark/abundance.txt \
	Embryo/abundance.txt \
	FlushBud/abundance.txt \
	MatureNeedle/abundance.txt \
	Megagametophyte/abundance.txt \
	SeedGermination/abundance.txt \
	Xylem/abundance.txt \
	YoungBuds/abundance.txt

salmon.html: \
	tissues.tsv \
	Bark_quant/quant.sf \
	Embryo_quant/quant.sf \
	FlushBud_quant/quant.sf \
	MatureNeedle_quant/quant.sf \
	Megagametophyte_quant/quant.sf \
	SeedGermination_quant/quant.sf \
	Xylem_quant/quant.sf \
	YoungBuds_quant/quant.sf

# Data

FlushBud.bam:
	ln -sf /projects/analysis/analysis7/IX0071/D0U8CACXX_4/no_ref/A11179/bwa/D0U8CACXX_4_TTAGGC.bam $@

Bark.bam:
	ln -sf /projects/analysis/analysis7/IX0071/D0U8CACXX_4/no_ref/A11180/bwa/D0U8CACXX_4_TGACCA.bam $@

Xylem.bam:
	ln -sf /projects/analysis/analysis7/IX0275/D0U8CACXX_6/no_ref/A11181/bwa/D0U8CACXX_6_ACAGTG.bam $@

MatureNeedle.bam:
	ln -sf /projects/analysis/analysis7/IX0275/D0U8CACXX_6/no_ref/A11182/bwa/D0U8CACXX_6_GCCAAT.bam $@

YoungBuds.bam:
	ln -sf /projects/analysis/analysis7/IX0072/D0U8CACXX_5/no_ref/A11183/bwa/D0U8CACXX_5_CAGATC.bam $@

Megagametophyte.bam:
	ln -sf /projects/analysis/analysis8/IX0729/D12M7ACXX_7/no_ref/A17447/bwa/D12M7ACXX_7_TAGCTT.bam $@

Embryo.bam:
	ln -sf /projects/analysis/analysis8/IX0729/D12M7ACXX_7/no_ref/A17448/bwa/D12M7ACXX_7_GGCTAC.bam $@

SeedGermination.bam:
	ln -sf /projects/analysis/analysis7/IX0072/D0U8CACXX_5/no_ref/A11276/bwa/D0U8CACXX_5_ACTTGA.bam $@
