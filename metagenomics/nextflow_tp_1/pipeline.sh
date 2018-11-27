#!/bin/bash

# Trimming
vsearch --fastq_convert 20ng-20cycles-1_R1.fastq.gz --fastq_qminout 20 --fastqout 20ng-20cycles-1_R1.trimmed.fastq
vsearch --fastq_convert 20ng-20cycles-1_R2.fastq.gz --fastq_qminout 20 --fastqout 20ng-20cycles-1_R2.trimmed.fastq

# Merging
vsearch --fastq_mergepairs 20ng-20cycles-1_R1.trimmed.fastq --reverse 20ng-20cycles-1_R2.trimmed.fastq --fastq_minovlen 40 --fastq_maxdiffs 15 --fastaout 20ng-20cycles-1.merged.fasta --label_suffix ";sample=20ng-20cycles-1"
awk '{gsub (" ", "", $0); print}' 20ng-20cycles-1.merged.fasta > 20ng-20cycles-1.merged.nospaces.fasta

# Merge of Merged


# Dereplication & Singleton removal
vsearch --derep_fulllength 20ng-20cycles-1.merged.nospaces.fasta --strand both --minuniquesize 10 --output 20ng-20cycles-1.dereplicated.fasta --sizeout

# Chimera filtering
vsearch --uchime_denovo 20ng-20cycles-1.dereplicated.fasta --nonchimeras 20ng-20cycles-1.nonchimeras.fasta

# Clustering
vsearch --cluster_size 20ng-20cycles-1.nonchimeras.fasta --id 0.97 --otutabout 20ng-20cycles-1.clusters.table --consout 20ng-20cycles-1.clusters.fasta

# Mapping
vsearch --usearch_global 20ng-20cycles-1.clusters.fasta --db mock_16S_18S.fasta --id 0.97 --alnout 20ng-20cycles-1.mock_16S_18S.alignment --samout 20ng-20cycles-1.mock_16S_18S.sam

# Print out species
awk '{print $3;}' 20ng-20cycles-1.mock_16S_18S.sam



