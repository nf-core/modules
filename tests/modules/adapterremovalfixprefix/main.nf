#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADAPTERREMOVAL          } from '../../../modules/adapterremoval/main.nf'
include { ADAPTERREMOVALFIXPREFIX } from '../../../modules/adapterremovalfixprefix/main.nf'

workflow test_adapterremovalfixprefix {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    ADAPTERREMOVAL ( input, [] )
    ADAPTERREMOVALFIXPREFIX ( ADAPTERREMOVAL.out.collapsed )
}
