# Written by Shaun Jackman @sjackman

# Parameters
genome=pg29mt-scaffolds
transcripts=transcripts
threads=4

# For shell redirection <(gunzip -c reads.fq.gz)
SHELL=/bin/bash

all: salmon.html rna-editing.html $(genome).vcf.gz

install-deps:
	brew install bcftools bwa datamash kallisto salmon samtools vcftools

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:

# Create symlinks

$(genome).gff:
	ln -sf ../$@ .

# Extract GFF attributes

%.gff.ID-Name.tsv: %.gff
	grep '\tgene\t' $< |gsed -E 's/.*ID=([^;]*);Name=([^;]*).*/\1\t\2/' >$@

# Convert BAM to FASTQ

%_1.fq.gz: %.bam
	time samtools view -uf64 $< |samtools bam2fq - |gzip >$@

%_2.fq.gz: %.bam
	time samtools view -uf128 $< |samtools bam2fq - |gzip >$@

# BWA

# Index the reference
%.fa.bwt: %.fa
	bwa index $<

# Align the paired FASTQ query to the genome
$(genome).%.sam: %_1.fq.gz %_2.fq.gz $(genome).fa.bwt
	bwa mem -t$(threads) $(genome).fa $*_1.fq.gz $*_2.fq.gz >$@

# Align the paired FASTQ query to the transcripts
$(transcripts).%.sam: %_1.fq.gz %_2.fq.gz $(transcripts).fa.bwt
	bwa mem -t$(threads) $(transcripts).fa $*_1.fq.gz $*_2.fq.gz >$@

# samtools

# Sort a SAM file and create a BAM file
%.sort.bam: %.sam
	samtools sort -@$(threads) -m8G -T/var/tmp/sjackman/$@ -Obam -o$@ $<

# Index a BAM file
%.bam.bai: %.bam
	samtools index $<

# Call variants
%.vcf.gz: %.fa
	samtools mpileup -vf $^ >$@

# bcftools

# Compute stats of a VCF file
%.vcf.gz.stats: %.vcf.gz
	bcftools stats $< >$@

# Call genotypes
%.gt.vcf.gz: %.vcf.gz
	bcftools call -mvf GQ,GP -Oz $< >$@

# Filter high-quality biallelic SNV
%.qual.vcf.gz: %.vcf.gz
	bcftools filter -S . -e 'QUAL < 50 || FMT/GQ < 50' $< |bcftools view -a |bcftools view -m2 -M2 -v snps -Oz >$@

# Identify regions of coding genes
%.regions.tsv: %.ID-Name.tsv
	egrep -v 'orf|rrn|trn' $< |awk '{print $$1 "\t1\t999999999"}' >$@

# Filter coding genes
%.cds.vcf.gz: %.vcf.gz $(genome).gff.regions.tsv
	bcftools view -T $(genome).gff.regions.tsv -Oz $< >$@

# Filter C-to-U edits
%.ctou.vcf.gz: %.vcf.gz
	bcftools view -v snps -i 'REF="C" && ALT~"T"' -Oz $< >$@

# Count the number of C-to-U edits per gene
%.ctou.tsv: %.vcf.gz.tsv
	(printf "ID\tC>U\n"; grep '\tC\t.*[CT]/T' $< |datamash -g 1 count 1) >$@

# vcftools

# Convert VCF to TSV

%.vcf.gz.tsv: %.vcf.gz
	gunzip -c $< |vcf-to-tab >$@

# tabix

# Index a bgz compressed VCF file
%.vcf.gz.tbi: %.vcf.gz
	tabix -pvcf $<

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
	mogrify -units PixelsPerInch -density 300 $*_files/figure-html/*.png

# Dependencies

$(genome).vcf.gz: \
	$(genome).Bark.sort.bam \
	$(genome).Embryo.sort.bam \
	$(genome).FlushBud.sort.bam \
	$(genome).MatureNeedle.sort.bam \
	$(genome).Megagametophyte.sort.bam \
	$(genome).SeedGermination.sort.bam \
	$(genome).Xylem.sort.bam \
	$(genome).YoungBuds.sort.bam

$(transcripts).vcf.gz: \
	$(transcripts).Bark.sort.bam \
	$(transcripts).Embryo.sort.bam \
	$(transcripts).FlushBud.sort.bam \
	$(transcripts).MatureNeedle.sort.bam \
	$(transcripts).Megagametophyte.sort.bam \
	$(transcripts).SeedGermination.sort.bam \
	$(transcripts).Xylem.sort.bam \
	$(transcripts).YoungBuds.sort.bam

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

rna-editing.html: transcripts.gt.qual.ctou.tsv

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
