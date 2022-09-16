#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAPPA_EXAMINEASSIGN } from '../../../../modules/gappa/examineassign/main.nf'

workflow test_gappa_examineassign {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GAPPA_EXAMINEASSIGN ( input )
}
