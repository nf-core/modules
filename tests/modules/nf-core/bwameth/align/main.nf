#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMETH_INDEX } from '../../../../modules/bwameth/index/main.nf'
include { BWAMETH_ALIGN } from '../../../../modules/bwameth/align/main.nf'

//
// Test with single-end data
//
workflow test_bwameth_align_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWAMETH_INDEX ( fasta )
    BWAMETH_ALIGN ( input, BWAMETH_INDEX.out.index )
}

//
// Test with paired-end data
//
workflow test_bwameth_align_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_methylated_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWAMETH_INDEX ( fasta )
    BWAMETH_ALIGN ( input, BWAMETH_INDEX.out.index )
}
