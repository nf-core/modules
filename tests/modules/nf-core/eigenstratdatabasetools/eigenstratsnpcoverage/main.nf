#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE } from '../../../../../modules/nf-core/eigenstratdatabasetools/eigenstratsnpcoverage/main.nf'

workflow test_eigenstratdatabasetools_eigenstratsnpcoverage {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE ( input )
}
