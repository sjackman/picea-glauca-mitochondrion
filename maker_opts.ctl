#-----Genome (these are always required)
genome=pg29mt-scaffolds.fa #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=cds_aa.fa #protein sequence file in fasta format (i.e. from mutiple oransisms)

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=picea #select a model organism for RepBase masking in RepeatMasker
rmlib=rmlib.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/usr/local/opt/maker/libexec/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner

#-----Gene Prediction
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
#trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no

#-----External Application Behavior Options
cpus=4 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
est_forward=1 #map names and attributes forward from EST evidence, 1 = yes, 0 = no
single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
#single_length=50 #min length required for single exon ESTs if 'single_exon is enabled'
