#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTK_FASTK } from '../../../../../modules/nf-core/fastk/fastk/main.nf'

workflow test_fastk_fastk_single_end {

    input = [
        [ id:'test' , single_end: true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    FASTK_FASTK ( input )
}

workflow test_fastk_fastk_paired_end {

    input = [
        [ id:'test' , single_end: false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FASTK_FASTK ( input )
}
