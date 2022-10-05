#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PEAR } from '../../../../modules/nf-core/pear/main.nf'

workflow test_pear {

    input = [
        [ id:'test' ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    PEAR ( input )
}
