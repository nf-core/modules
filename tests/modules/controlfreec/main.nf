#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONTROLFREEC } from '../../../modules/controlfreec/main.nf'

workflow test_controlfreec {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    CONTROLFREEC ( input )
}
