#!/bin/bash

# Uncompress
# gunzip -vk *.fastq.gz

# Trimming
# alientrimmer -if EchA_R1.fastq -ir EchA_R2.fastq -q 20 -l 45 -c ../soft/AlienTrimmer_0.4.0/alienTrimmerPF8contaminants.fasta -of EchA_R1.trimmed.fastq -or EchA_R2.trimmed.fastq

# Assembly
megahit -1 EchA_R1.trimmed.fastq -2 EchA_R2.trimmed.fastq -m 0.5 -t 3 --k-list 39 -o EchA_megahit

# Gene prediction
prodigal -i EchA_megahit/final.contigs.fa -p meta -c -o EchA.prodigal.wholegenes.gbk -a EchA.prodigal.wholegenes.proteins.fa
awk '{gsub (" ", "", $0); print}' EchA.prodigal.wholegenes.proteins.fa > EchA.prodigal.wholegenes.proteins.nospaces.fa

# Clustering
cd-hit -i EchA.prodigal.wholegenes.proteins.nospaces.fa -o EchA.cd-hit -M 4000 -T 3 -aS 0.9 -c 0.95 -d 0

# Blastp
makeblastdb -in resfinder.faa -dbtype prot
blastp -query EchA.cd-hit -out EchA.blastp -outfmt "6 qseqid sseqid qcovs" -task blastp -db resfinder.faa -evalue 1e-5 -num_threads 3 -qcov_hsp_perc 80 -max_target_seqs 1

 







