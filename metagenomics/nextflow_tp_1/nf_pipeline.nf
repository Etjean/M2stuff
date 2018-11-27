// Pipeline en Nextflow


params.in = "$HOME/M2stuff/metagenomics/tp_1/my_dvpt2/20ng-20cycles-1_R1.fastq.gz"

reads = file(params.in)


process trimReads {

    input:
    file input from reads

    output:
    file output

    """
    vsearch --fastq_convert $input --fastq_qminout 20 --fastqout output
    """

}



// process reverse {

//     input:
//     file x from records
    
//     output:
//     stdout result

//     """
//     cat $x | rev
//     """
// }



println 'c\'est fini'