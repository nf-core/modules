#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_INDEX }                     from '../../../../modules/biscuit/index/main.nf'
include { BISCUIT_ALIGN as BISCUIT_ALIGN_SE } from '../../../../modules/biscuit/align/main.nf'
include { BISCUIT_ALIGN as BISCUIT_ALIGN_PE } from '../../../../modules/biscuit/align/main.nf'


// Single-end test
workflow test_biscuit_align_single {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX ( fasta )
    BISCUIT_ALIGN_SE (input, BISCUIT_INDEX.out.index )
}

// paired-end test
workflow test_biscuit_align_paired {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_methylated_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX ( fasta )
    BISCUIT_ALIGN_SE (input, BISCUIT_INDEX.out.index )
}
