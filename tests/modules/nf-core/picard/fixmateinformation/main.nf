#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_FIXMATEINFORMATION } from '../../../../modules/picard/fixmateinformation/main.nf'

workflow test_picard_fixmateinformation {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    PICARD_FIXMATEINFORMATION ( input )
}
