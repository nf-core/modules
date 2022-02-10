#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_BISCUITBLASTER } from '../../../../modules/biscuit/biscuitblaster/main.nf'

workflow test_biscuit_biscuitblaster {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    BISCUIT_BISCUITBLASTER ( input )
}
