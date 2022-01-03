#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_ALIGN } from '../../../../modules/biscuit/align/main.nf'

workflow test_biscuit_align {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    BISCUIT_ALIGN ( input )
}
