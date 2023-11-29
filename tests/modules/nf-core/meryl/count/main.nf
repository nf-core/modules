#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERYL_COUNT } from '../../../../../modules/nf-core/meryl/count/main.nf'

workflow test_meryl_count_single_end {

    input = [
        [ id:'test' , single_end: true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MERYL_COUNT ( input )
}

workflow test_meryl_count_paired_end {

    input = [
        [ id:'test' , single_end: false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    MERYL_COUNT ( input )
}
