#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERYL_COUNT    } from '../../../../../modules/nf-core/meryl/count/main.nf'
include { MERYL_UNIONSUM } from '../../../../../modules/nf-core/meryl/unionsum/main.nf'

workflow test_meryl_unionsum_single_end {

    input = [
        [ id:'test', single_end: true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MERYL_COUNT ( input )
    MERYL_UNIONSUM ( MERYL_COUNT.out.meryl_db )
}

workflow test_meryl_unionsum_paired_end {

    input = [
        [ id:'test', single_end: false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    MERYL_COUNT ( input )
    MERYL_UNIONSUM ( MERYL_COUNT.out.meryl_db )
}
