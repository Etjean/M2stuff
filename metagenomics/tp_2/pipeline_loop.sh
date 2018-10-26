#!/bin/bash

dirfq=fastq/
dirdb=my_dvpt/
dir=results/

# Uncompress
# gunzip -vk *.fastq.gz

# Trimming
# alientrimmer -if EchA_R1.fastq -ir EchA_R2.fastq -q 20 -l 45 -c ../soft/AlienTrimmer_0.4.0/alienTrimmerPF8contaminants.fasta -of EchA_R1.trimmed.fastq -or EchA_R2.trimmed.fastq


for f in $dirfq/*_R1.trimmed.fastq ; do
    # Get basename
    fbn=$(basename $f _R1.trimmed.fastq)

    # # Assembly
    # megahit -1 $dirfq$fbn'_R1.trimmed.fastq' -2 $dirfq$fbn'_R2.trimmed.fastq' -m 0.6 -t 3 --k-list 39 -o $dir$fbn'_megahit'

    # # Gene prediction
    # prodigal -i $dir$fbn'_megahit/final.contigs.fa' -p meta -c -a $dir$fbn'.prodigal.wholegenes.proteins.fa'
    # awk '{gsub (" ", "", $0); print}' $dir$fbn'.prodigal.wholegenes.proteins.fa' > $dir$fbn'.prodigal.wholegenes.proteins.nospaces.fa'

    # Merging all fasta proteins files
    cat $dir$fbn'.prodigal.wholegenes.proteins.nospaces.fa' >> $dir'all.prodigal.wholegenes.proteins.nospaces.fa'

done




# Clustering
cd-hit -i $dir'all.prodigal.wholegenes.proteins.nospaces.fa' -o $dir'all.cd-hit' -M 4000 -T 3 -aS 0.9 -c 0.95 -d 0

# Blastp
# makeblastdb -in $dirdb'resfinder.faa' -dbtype prot
blastp -query $dir'all.cd-hit' -out $dir'all.blastp' -outfmt "6 qseqid sseqid qcovs" -task blastp -db $dirdb'resfinder.faa' -evalue 1e-5 -num_threads 3 -qcov_hsp_perc 80 -max_target_seqs 1