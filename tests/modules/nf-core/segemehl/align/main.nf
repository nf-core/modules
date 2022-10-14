#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEGEMEHL_ALIGN } from '../../../../../modules/nf-core/segemehl/align/main.nf'

workflow test_segemehl_align_pe {

    input = [
        [ id:'test', single_end:false ],

        [
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    index = file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/segemehl/genome.idx")

    SEGEMEHL_ALIGN ( input, fasta, index )
}

workflow test_segemehl_align_se {

    input = [
        [ id:'test', single_end:true ],

        [ file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true) ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    index = file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/segemehl/genome.idx")

    SEGEMEHL_ALIGN ( input, fasta, index )
}
