# Annotate and visualize the white spruce mitochondrial genome
# Copyright 2015 Shaun Jackman

# Name of the assembly
name=pg29mt-scaffolds

# Number of threads
t=4

# Parallel compression with pigz.
gzip=pigz -p$t

# Green plant mitochondria
edirect_query='Viridiplantae[Organism] mitochondrion[Title] (complete genome[Title] OR complete sequence[Title])'

all: $(name).gff $(name).gbk $(name).tbl $(name).sqn \
	$(name).maker.evidence.gff $(name).maker.repeat.gff \
	$(name).maker.gff.gene $(name).prokka.gff.gene $(name).gff.gene $(name).tbl.gene \
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

# Align the scaffolds to Cycas taitungensis NC_010303.1
NC_010303.1.%.blastn: %.fa NC_010303.1.fa
	blastn -subject NC_010303.1.fa -query $< -out $@

# Align proteins to the NR database and output TSV.
%.blastp.tsv: %.fa
	(printf 'qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tstaxid\tstitle\n'; \
	blastp -remote -db nr -max_target_seqs 5 -outfmt '6 std qlen slen qcovs staxid stitle' -query $<) >$@

# Align proteins to the NR database and output PAF.
%.blastp.paf.gz: %.fa
	blastp -remote -db nr -max_target_seqs 5 -outfmt '6 qaccver qlen qstart qend sstrand saccver slen sstart send nident length evalue staxid stitle' -query $< \
	| awk -F'\t' -vOFS='\t' '{ $$5 = "+"; $$13 = "tx:i:" $$13; $$14 = "ti:z:" $$14; print }' \
	| $(gzip) >$@

# BWA
#-------------------------------------------------------------------------------

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align a FASTA file to the draft genome using BWA-MEM.
$(ref).bwa.%.sam.gz: %.fa $(ref).fa.bwt
	bwa mem $(ref).fa $< | $(gzip) >$@

# minimap2
#-------------------------------------------------------------------------------

# Align a FASTA file to the draft genome using minimap2.
$(ref).minimap2.%.sam.gz: $(ref).fa %.fa
	minimap2 -a -xsplice $^ | $(gzip) >$@

# Copy local data
#-------------------------------------------------------------------------------

#PICEAGLAUCA_rpt2.0.fa: /genesis/extscratch/seqdev/PG/data/PICEAGLAUCA_rpt2.0
#	cp -a $< $@

# NCBI
#-------------------------------------------------------------------------------

# Download the Picea glauca plastid FASTA.
pglaucacp.fa:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id=NC_028594.1&retmode=text&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Download the Picea glauca plastid CDS FASTA.
pglaucacp.cds.orig.fa:
	curl -o $@ 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id=NC_028594.1&retmode=text&db=nucleotide&rettype=fasta_cds_na'

# Download the Picea glauca plastid GFF.
pglaucacp.gff:
	curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?id=NC_028594.1&db=nuccore&report=gff3' >$@

# Download the Picea glauca mitochondrion FASTA.
pg29mt-scaffolds.orig.fa:
	curl ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LK/AM/LKAM01/LKAM01.1.fsa_nt.gz | gunzip -c | seqtk seq >$@

# Rename the sequences.
pg29mt-scaffolds.fa: pg29mt-scaffolds.orig.fa
	gsed -E 's/^>(.*gb[|]([^|]*).*)[|]/>\2 \1/' $< >$@

# Download the Picea glauca mitochondrion GFF.
LKAM01.1.gff:
	(echo '##gff-version 3'; \
	for i in $$(seq -w 1 36); do \
		curl -s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=LKAM010000$$i.1"; \
	done) | sed '/^$$/d' >$@

# Download the Pinus strobus CDS FASTA.
pstrobusmt.cds.orig.fa:
	esearch -db nuccore -query 'Pinus strobus[Organism] AND mitochondrion[filter]' \
		| efetch -format fasta_cds_na \
		| seqmagick convert --pattern-include 'AJP335' --line-wrap=0 - $@

# Download the Pinus strobus protein FASTA.
pstrobusmt.aa.orig.fa:
	esearch -db nuccore -query 'Pinus strobus[Organism] AND mitochondrion[filter]' \
		| efetch -format fasta_cds_aa \
		| seqmagick convert --pattern-include 'AJP335' --line-wrap=0 - $@

# Rename the coding sequences.
%.cds.fa: %.cds.orig.fa
	gsed -E \
		-e 's/>(.*gene=([^]]*)].*protein_id=([^]]*)].*)/>\2_\3 \1/' -e '/gene=.*protein_id=/n' \
		-e 's/>(.*gene=([^]]*)].*locus_tag=([^]]*)].*)/>\2_\3 \1/' -e '/gene=.*locus_tag=/n' \
		-e 's/>(.*protein_id=([^]]*)].*)/>\2 \1/' -e '/protein_id=/n' \
		-e 's/>(.*locus_tag=([^]]*)].*)/>\2 \1/' $< \
	| seqmagick convert --deduplicate-sequences --line-wrap=0 - $@
	seqmagick convert --sort=name-asc --line-wrap=0 $@ $@

# Rename the protein sequences.
%.aa.fa: %.aa.orig.fa
	gsed -E \
		-e 's/>(.*gene=([^]]*)].*protein_id=([^]]*)].*)/>\2_\3 \1/' -e '/gene=.*protein_id=/n' \
		-e 's/>(.*gene=([^]]*)].*locus_tag=([^]]*)].*)/>\2_\3 \1/' -e '/gene=.*locus_tag=/n' \
		-e 's/>(.*protein_id=([^]]*)].*)/>\2 \1/' -e '/protein_id=/n' \
		-e 's/>(.*locus_tag=([^]]*)].*)/>\2 \1/' $< \
	| seqmagick convert --deduplicate-sequences --line-wrap=0 - $@
	seqmagick convert --sort=name-asc --line-wrap=0 $@ $@

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
		--pattern-replace 'gene=RNA_pol' 'gene=rpo' \
		--pattern-replace 'gene=yejU' 'gene=ccmC' \
		--pattern-replace 'gene=yejV' 'gene=ccmB' \
		--pattern-replace 'gene=18S rRNA' 'gene=40' \
		$< - \
	|gsed -E '/protein=[^]]*intron[^]]*ORF/s/gene=/gene=ymf/; \
		s/^>(.*gene=([^]]*).*)$$/>\2|\1/' \
	|seqmagick -q convert --pattern-exclude '^lcl' --deduplicate-taxa - $@

# Extract accession numbers from the FASTA file
%.id: %.orig.fa
	sed '/^>/!d;s/.*lcl|//;s/_prot_.*//' cds_aa.orig.fa |uniq |sort -u >$@

# Fetch the records
%.docsum.xml: %.id
	esearch -db nuccore -query "`<$<`" |efetch -format docsum >$@

# Convert XML to TSV
%.docsum.tsv: %.docsum.xml
	(printf "Caption\tTaxId\tOrganism\tTitle\n"; \
		xtract -pattern DocumentSummary -element Caption,TaxId,Organism,Title <$<) >$@

# Fetch a taxonomy tree from a list of TaxID
%.taxid.xml: %.docsum.tsv
	awk 'NR>1 {print $$2}' $< |sort -u |epost -db taxonomy |efetch -format xml >$@

# Convert NCBI XML taxonomy to Newick
%.taxid.tree: %.taxid.xml
	xsltproc bin/taxon2newick.xsl $< >$@

# Cycas taitungensis
NC_010303.1.json: %.json:
	bionode-ncbi search nuccore $* >$@

# Extract a list of a gene names
NC_010303.1.gff.gene: NC_010303.1.gff
	gsed -nE 's/ND/nad/;/\tgene\t/s/.*ID=([^;]*).*/\1/p' $< |sort -u >$@

# Extract a list of tRNA gene names
NC_010303.1.gff.tRNA: NC_010303.1.gff
	gsed -nE '/\ttRNA\t/!d;s/.*ID=([^;]*).*codon_recognized=([^;]*).*/\1-\2/p' $< \
		|bioawk -F- '{print $$1 "-" revcomp($$2)}' |tr T U \
		|sort -u >$@

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

# Samtools
#-------------------------------------------------------------------------------

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file by target position.
%.sort.bam: %.sam.gz
	samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index -@$t $<

# Prodigal
#-------------------------------------------------------------------------------

# Annotate genes using Prodigal
%.prodigal.gff: %.fa
	prodigal -c -m -g 1 -p single -f gff -a $*.prodigal.faa -d $*.prodigal.ffn -s $*.prodigal.tsv -i $< -o $@

# RepeatModeler
#-------------------------------------------------------------------------------

%.nin: %.fa
	BuildDatabase -name $* -engine ncbi $<

%.RepeatModeler.fa: %.nin
	RepeatModeler -database $*
	cp -a RM_*/consensi.fa.classified $@

# ARAGORN
#-------------------------------------------------------------------------------

# Annotate tRNA using ARAGORN and output TSV
%.aragorn.tsv: %.fa
	aragorn -gcstd -l -w -o $@ $<

# Annotate tRNA using ARAGORN and output text
%.aragorn.txt: %.fa
	aragorn -gcstd -l -o $@ $<

# Convert ARAGORN output to GFF
%.aragorn.gff: %.aragorn.tsv
	bin/aragorn_out_to_gff3.py --full <$< |gt gff3 -sort |bin/gt-bequeath Name |grep -v trnX >$@

# Annotate tRNA using tRNAscan-SE
%.trnascan.orig.tsv: %.fa
	tRNAscan-SE -O -o $@ -f $*.trnascan.txt $<

# Barrnap
#-------------------------------------------------------------------------------

# Identify ribosomal RNA (rRNA) genes.
%.barrnap.orig.gff: %.fa
	barrnap --kingdom bac --threads $t $< >$@

# Rename the rRNA genes.
%.barrnap.gff: %.barrnap.orig.gff
	sed -e 's/16S_rRNA/rrn16/' -e 's/23S_rRNA/rrn23/' $< >$@

# RNAmmer
#-------------------------------------------------------------------------------

# Annotate rRNA using RNAmmer
%.rnammer.gff2: %.fa
	mkdir -p rnammer
	rnammer -S bac -gff $@ -xml rnammer/$*.xml -f rnammer/$*.fa -h rnammer/$*.hmm $<

# Convert GFF2 to GFF3
%.rnammer.gff: %.rnammer.gff2
	bin/convert_RNAmmer_to_gff3.pl --input=$< \
		|sed -E -e 's/ID=([^s]*)s_rRNA_([0-9]*)/Name=rrn\1;&/' \
			-e 's/rrn16/rrn18/g;s/rrn23/rrn26/g' >$@

# Infernal and RFAM
#-------------------------------------------------------------------------------
# See https://rfam.readthedocs.io/en/latest/genome-annotation.html

# Download RFAM clans.
Rfam.clanin:
	curl -o $@ ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.0/Rfam.clanin

# Download RFAM.
Rfam.cm:
	curl ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.0/Rfam.cm.gz | gunzip -c >$@

# Download RFAM motif Domain-V (RM00007).
RM00007.cm:
	curl -o $@ http://rfam.xfam.org/motif/RM00007/cm

# Add the accession number to domain V.
domainV.cm: RM00007.cm
	gsed 's/^NAME.*$$/&\nACC      RM00007/' $< >$@

# Download RFAM covariance models.
RF%.cm:
	curl -o $@ http://rfam.xfam.org/family/RF$*/cm

# Download RFAM family Intron_gpII (RF00029), which includes domains V and VI.
# Decrease the bit score threshold from 40 to 20.
groupII.cm: RF00029.cm
	sed 's/^GA.*/GA       20.00/' $< >$@

# Download RFAM clan group-II-D1D4 (CL00102).
# Decrease the bit score threshold from to 10.
d1d4.cm: RF01998.cm RF01999.cm RF02001.cm RF02003.cm RF02004.cm RF02005.cm RF02012.cm
	sed 's/^GA.*/GA       10.00/' $^ >$@

# Download a covariance model of bacterial 5S, 18S, and 26S rRNA.
rRNA.cm: RF00001.cm RF00177.cm RF02541.cm
	cat $^ >$@

# Download a covariance model of tRNA.
tRNA.cm: RF00005.cm
	cat $^ >$@

# Download a covariance model of RFAM families relevant to mitochondria.
mt.cm: rRNA.cm tRNA.cm domainV.cm groupII.cm d1d4.cm
	cat $^ >$@

# Create a clan file for the mitochondrial covariance model.
mt.clanin: Rfam.clanin
	grep CL00102 $< >$@

# Compress RFAM using Infernal.
%.cm.i1m: %.cm
	cmpress $<

# Identify RFAM families for mitochondria using Infernal.
%.mt.infernal: %.fa mt.clanin mt.cm.i1m
	cmscan --cpu=$t --mid --cut_ga --fmt 2 --tblout $*.mt.infernal --clanin mt.clanin mt.cm $< >$*.mt.txt

# Identify RFAM families using Infernal.
%.rfam.infernal: %.fa Rfam.clanin Rfam.cm.i1m
	cmscan --cpu=$t --mid --cut_ga --fmt 2 --tblout $*.rfam.infernal --clanin Rfam.clanin Rfam.cm $< >$*.rfam.txt

# Identify group II domain V motif using Infernal.
%.domainV.infernal: %.fa domainV.cm.i1m
	cmscan --cpu=$t --mid --cut_ga --fmt 2 --tblout $*.domainV.infernal domainV.cm $< >$*.domainV.txt

# Identify group II domains V and VI family using Infernal.
%.groupII.infernal: %.fa groupII.cm.i1m
	cmscan --cpu=$t --mid --incT=20 --fmt 2 --tblout $@ groupII.cm $< >$*.groupII.txt

# Identify RFAM clan group-II-D1D4 (CL00102) using Infernal.
%.d1d4.infernal: %.fa d1d4.cm.i1m
	cmscan --cpu=$t --nohmm --incT=10 --fmt 2 --tblout $@ d1d4.cm $< >$*.d1d4.txt

# Convert Infernal TBL format to GFF.
%.gff: %.infernal
	awk -vOFS='\t' 'BEGIN { print "##gff-version 3" } /^#/ { next } { a = $$10; b = $$11; } $$10 >= $$11 { a = $$11; b = $$10 } { print $$4, "Infernal", "match", a, b, $$17, $$12, ".", "Name=" $$2 ";Target=" $$3 " " $$8 " " $$9 ";e=" $$18 }' $< \
	| gt gff3 -sort >$@

# MAKER
#-------------------------------------------------------------------------------

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
	| gsed -E \
		-e '/\tintron\t/d' \
		-e 's/Name=trnascan-[^-]*-noncoding-([^-]*)-gene/Name=trn\1/g' \
		-e '/\trRNA\t/s/ID=([^;]*)s_rRNA/Name=rrn\1;&/g' \
	| gt gff3 -addintrons -sort - >$@

# Add the rRNA annotations to the GFF file
$(name).maker.gff: $(name).rnammer.gff

# Add the tRNA annotations to the GFF file
$(name).maker.gff: $(name).aragorn.gff

# Prokka
#-------------------------------------------------------------------------------

# Convert the FASTA file to the Prokka FASTA format
cds_aa.prokka.fa: %.prokka.fa: %.fa
	sed -E 's/^>([^ ]*) .*gene=([^]]*).*protein=([^]]*).*$$/>\1 ~~~\2~~~\3/; s/^-//' $< >$@

# Annotate genes using Prokka
prokka/%.gff: %.fa cds_aa.prokka.fa
	prokka --kingdom bac --gcode 1 --addgenes --proteins cds_aa.prokka.fa --rnammer \
		--cpus $t \
		--genus Picea --species 'glauca mitochondrion' \
		--locustag PG29MT \
		--force --outdir prokka --prefix $* \
		$<

# Remove the FASTA section from the Prokka GFF file
%.prokka.gff: prokka/%.gff
	gsed -E -e '/^##FASTA/,$$d' \
		-e 's/=atpD/=atpB/g' \
		-e 's/=ltrA/=ymfltrA/g' \
		-e 's/=rplB/=rpl2/g' \
		-e 's/=smc/=dpo/g' \
		-e '/\tgene\t/{/gene=/!s/ID=[^_]*_([0-9]*)/Name=orf\1;&/;}' \
		-e '/\tCDS\t/{/gene=/!s/ID=[^_]*_([0-9]*)/Name=orf\1;&/; s/CDS/mRNA/;p; s/mRNA/CDS/;s/Parent=[^;]*;//;s/ID=/Parent=/;}' \
		$< >$@

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
	(cat $< && sed -En -e '/^FEATURES/,$$ { s/Name="([^|]*).*"/gene="\1"/; p; }' $*.gb) >$@

# Merge MAKER and Prokka annotations using bedtools
%.maker.prokka.gff: %.prokka.gff %.maker.gff
	bedtools intersect -v -header -a $< -b $*.maker.gff \
		|sed '/tRNA-???/{N;d;}' \
		|gt gff3 -sort $*.maker.gff - >$@

# Merge manual, MAKER and Prokka annotations using bedtools
%.orig.gff: %.maker.prokka.gff %.manual.gff
	bedtools intersect -v -header -a $< -b $*.manual.gff \
		|sed '/tRNA-???/{N;d;}' \
		|gt gff3 -sort $*.manual.gff - >$@

# Remove contamination from the GFF file
%.gff: %.orig.gff
	sed -E '/^(33|36)[[:blank:]]/d' $< |uniq >$@

# LKAM01.2 Picea glauca mitochondrion
#-------------------------------------------------------------------------------

# Merge the previous annotation LKAM01.1 with the current manual annotation.
LKAM01.2.gff: LKAM01.1.gff pg29mt-scaffolds.gff
	bedtools intersect -v -header -a $< -b pg29mt-scaffolds.manual.gff \
		| gt gff3 -sort pg29mt-scaffolds.manual.gff - >$@

# OrganellarGenomeDRAW
#-------------------------------------------------------------------------------

# Draw a linear genome
%.gbk.png: %.gbk
	drawgenemap --density 150 --format png --infile $< --outfile $<

# Draw a circular genome
%.gbk.circular.png: %.gbk
	drawgenemap --format png --infile $< --outfile $<.circular \
		--gc --force_circular --density 126
	mogrify -units PixelsPerInch -density 300 $@

# Report the annotated genes

# Convert a GFF file to a TSV file
%.gff.tsv: %.gff
	(printf "SeqID\tSource\tType\tStart\tEnd\tStrand\tName\n"; \
		egrep '\t(mRNA|rRNA|tRNA)\t' $< \
		|cut -sf 1,2,3,4,5,7,9 \
		|gsed 's/\t[^\t]*Name=/\t/;s/;.*//;s/|.*//;s/[-_].$$//;s/Prodigal:[0-9.]*/prokka/' \
		) >$@

# Extract the names of genes from a GFF file
%.gff.gene: %.gff
	bin/gff-gene-name $< >$@

# Determine unique gene names
%.gff.uniq.gene: %.gff.gene
	sed '/orf/d;s/|.*//;s/[-_][0-9a-z]$$//' $< |sort -u >$@

# Extract DNA sequences of GFF gene features from a FASTA file
%.gff.gene.fa: %.gff %.fa
	gt extractfeat -type gene -coords -matchdescstart -retainids -seqid -seqfile $*.fa $< >$@

# Extract DNA sequences of GFF CDS features from a FASTA file
%.gff.cds.fa: %.gff %.fa
	gsed -E 's/Name=([^;]*)/Target=\1 0 0/' $< \
	| gt extractfeat -type CDS -join -coords -target -matchdescstart -retainids -seqid -seqfile $*.fa - \
	| gsed -E 's/>(.*target IDs ([^]|]*).*)/>\2_\1/' >$@

# Extract the coding sequences of known genes.
%.cds.known.fa: %.cds.fa
	seqmagick convert --sort=name-asc --line-wrap=0 \
		 --pattern-include='atp1|atp4|atp6|atp8|atp9|ccmB|ccmC|ccmFC|ccmFN|cob|cox1|cox2|cox3|dpo|matR|mttB|nad1|nad2|nad3|nad4L|nad4|nad5|nad6|nad7|nad9|rpl10|rpl16|rpl2|rpl5|rpo|rps10|rps11|rps12|rps13|rps14|rps19|rps1|rps2|rps3|rps4|rps7|sdh3|sdh4' \
		 $< - \
	| sed -e 's/_part1//g' -e '/_part[2-9]/d' \
	| seqtk seq >$@

# Extract aa sequences of GFF CDS features from a FASTA file
%.gff.aa.fa: %.gff %.fa
	gsed -E 's/^.*Parent=([^;]*).*Name=([^;]*).*/&;Target=\2_\1 0 0/' $< \
	| gt extractfeat -type CDS -join -translate -coords -target -matchdescstart -retainids -seqid -seqfile $*.fa >$@

# Translate protein sequences of GFF CDS features from a FASTA file
%.aa.fa: %.fa
	gt seqtranslate -reverse no -fastawidth 0 $< \
	|sed -n '/ (1+)$$/{s/ (1+)$$//;p;n;p;n;n;n;n;}' >$@

# Extract sequences of GFF intron features
%.gff.intron.fa: %.gff %.fa
	gsed -E 's/Name=([^;]*)/Target=\1 0 0/' $< \
	| gt extractfeat -type intron -coords -target -matchdescstart -retainids -seqid -seqfile $*.fa - \
	| gsed -E 's/>(.*target IDs ([^]|]*).*)/>\2 \1/' >$@

# Extract sequences of GFF intron features plus flanking sequence
%.gff.intron.flank100.fa: %.gff %.fa.fai
	(awk '$$3 == "intron"' $<; \
		awk '$$3 == "intron"' $< |bedtools flank -b 100 -i stdin -g $*.fa.fai) \
	|sort -k1,1n -k4,4n - \
	|bedtools merge -i stdin \
	|bedtools getfasta -bed stdin -fi $*.fa -fo $@

# Extract sequences of GFF rRNA features
%.gff.rRNA.fa: %.gff %.fa
	gt extractfeat -type rRNA -coords -matchdescstart -retainids -seqid -seqfile $*.fa $< >$@

# Extract sequences of GFF tRNA features
%.gff.tRNA.fa: %.gff %.fa
	gt extractfeat -type tRNA -coords -matchdescstart -retainids -seqid -seqfile $*.fa $< >$@

# UniqTag

# Generate UniqTag from DNA or amino acid sequence
%.uniqtag: %.fa
	uniqtag $< >$@

# GenBank

# Split the GenBank file into one sequence per file
gbk/%.00.gbk: %.gbk
	gcsplit -sz -f $*. --suppress-matched $< '/\/\//' '{*}'
	mkdir -p $(@D)
	for i in $*.[0-9][0-9]; do mv $$i gbk/$$i.gbk; done

# Combine the OGDraw images into a single image
%.gbk.montage.png: \
		gbk/%.00.gbk.png \
		gbk/%.01.gbk.png \
		gbk/%.02.gbk.png \
		gbk/%.03.gbk.png \
		gbk/%.04.gbk.png \
		gbk/%.05.gbk.png \
		gbk/%.06.gbk.png \
		gbk/%.07.gbk.png \
		gbk/%.08.gbk.png \
		gbk/%.09.gbk.png \
		gbk/%.10.gbk.png \
		gbk/%.11.gbk.png \
		gbk/%.15.gbk.png \
		gbk/%.16.gbk.png \
		gbk/%.17.gbk.png \
		gbk/%.20.gbk.png \
		gbk/%.23.gbk.png
	montage -tile 3 -geometry +0+0 -units PixelsPerInch -density 1200 $^ gbk/$*.00.gbk_legend.png $@

# GenomeTools sketch
%.gff.png: %.gff
	gt sketch $@ $<

# Convert GFF to GTF
%.gtf: %.gff
	gt -q gff3_to_gtf $< >$@

# Convert GFF to TBL
%.tbl: %.gff %.product.tsv %.gff.aa.fa
	bin/gff3-to-tbl --centre=BCGSC --locustag=PG29MT $^ >$@

# Extract the names of genes from a TBL file
%.tbl.gene: %.tbl
	awk '$$1 == "gene" {print $$2}' $< >$@

# Extract the names of genes from MFannot annotation
%.mfannot.gene: %.mfannot
	gsed '/List of genes added/,/^;;$$/!d;s/^;; *//;/^-/d;/^$$/d;s/  *$$//;s/   */\n/g;s/ .*//' $< \
		|gsed '/orf/d;s/(/-/;s/-.*/\U&/;s/)//' >$@

# Remove contamination and add structured comments to the FASTA file
%.fsa: %.fa
	seqtk seq $< \
		|sed -E -e '/^>(33|36)$$/{N;d;}' \
			-e 's/^>.*/& [organism=Picea glauca] [isolate=PG29] [location=mitochondrion] [completeness=draft] [topology=linear] [gcode=1]/' \
			>$@

# tbl2asn
#-------------------------------------------------------------------------------

# Extract the IDs of protein sequences with unusual start codons.
%.alt: %.gff
	grep 'note=.*start codon' $< | egrep -o 'mRNA[0-9]+' >$@

# Extract protein sequences with unusual start codons.
%.pep: %.gff.aa.fa %.alt
	seqtk subseq $^ \
	| sed -E -e 's/^>(.*(gene[0-9]*).*)/>gnl|BCGSC|PG29MT_\2 \1/' \
		-e 's/^T/M/' \
		-e 's/^A/V/' >$@

# Convert TBL to GBK and SQN.
%.gbf %.sqn: %.fsa %.sbt %.tbl %.cmt %.pep
	tbl2asn -a s -i $< -t $*.sbt -w $*.cmt -Z $*.discrep -Vbv
	gsed -i 's/DEFINITION  Picea glauca/& mitochondrion draft genome/' $*.gbf

# Render HTML from RMarkdown.
%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<", output_format = "html_document")'
	mogrify -units PixelsPerInch -density 300 $*_files/figure-html/*.png

# Dependencies

genes.html: pg29mt-scaffolds.gff

repeats.html: pg29mt-scaffolds.maker.repeat.gff
