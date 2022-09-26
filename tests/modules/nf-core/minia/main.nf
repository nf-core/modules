#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIA } from '../../../modules/minia/main.nf'

workflow test_minia {
    input = [ [ id:'test' ], // meta map
              [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)]
            ]

    MINIA ( input )
}
