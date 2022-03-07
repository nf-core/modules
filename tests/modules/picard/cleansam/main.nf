#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_CLEANSAM } from '../../../../modules/picard/cleansam/main.nf'

workflow test_picard_cleansam {
    
    input = [ 
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) 
    ]

    PICARD_CLEANSAM ( input )
}
