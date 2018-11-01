# #!/bin/bash

# Input directories
dirfq='fastq/'
dirdb='databases/'
# Result directory
dir='results/'

# Trimming
for fq in $dirfq/* ; do
    # Get basename
    fbname=$(basename $fq .fastq.gz)
    # Trimming
    vsearch --fastq_convert $fq --fastq_qminout 20 --fastqout $dir$fbname".trimmed.fastq"
done


# Merging
for f in results/*_R1.trimmed.fastq ; do 
    # Get basename
    fbname=$(basename $f _R1.trimmed.fastq)
    # Merging
    vsearch --fastq_mergepairs $dir$fbname'_R1.trimmed.fastq' --reverse $dir$fbname'_R2.trimmed.fastq' --fastq_minovlen 40 --fastq_maxdiffs 15 --fastaout $dir$fbname'.merged.fasta' --label_suffix ';sample='$fbname
    awk '{gsub (" ", "", $0); print}' $dir$fbname'.merged.fasta' > $dir$fbname'.merged.nospaces.fasta'
    # Merging all merged files
    cat $dir$fbname'.merged.nospaces.fasta' >> $dir'all_merged.fasta'
done


# Dereplication & Singleton removal
vsearch --derep_fulllength $dir'all_merged.fasta' --strand both --minuniquesize 10 --output $dir'all_merged.dereplicated.fasta' #--sizeout

# Chimera filtering
vsearch --uchime_denovo $dir'all_merged.dereplicated.fasta' --nonchimeras $dir'all_merged.nonchimeras.fasta'

Clustering
vsearch --cluster_size $dir'all_merged.nonchimeras.fasta' --id 0.97 --consout $dir'all_merged.clusters.fasta' --otutabout $dir'all_merged.clusters.table'

# Mapping
vsearch --usearch_global $dir'all_merged.clusters.fasta' --db $dirdb'mock_16S_18S.fasta' --id 0.97 --samout $dir'all_merged.mock_16S_18S.sam' --alnout $dir'all_merged.mock_16S_18S.alignment' --userfields 'query+target' --userout $dir'all_merged.mock_16S_18S.annot_table'

# Print out species
echo ''
echo 'OTUS & SPECIES FOUND :'
cat $dir'all_merged.mock_16S_18S.annot_table'