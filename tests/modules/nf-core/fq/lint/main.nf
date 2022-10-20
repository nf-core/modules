#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FQ_LINT } from '../../../../../modules/nf-core/fq/lint/main.nf'

workflow test_fq_lint_success {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FQ_LINT ( input )
}

workflow test_fq_lint_fail {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FQ_LINT ( input )
}
