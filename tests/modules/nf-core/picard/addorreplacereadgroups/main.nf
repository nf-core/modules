#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_ADDORREPLACEREADGROUPS } from '../../../../modules/picard/addorreplacereadgroups/main.nf'

workflow test_picard_addorreplacereadgroups {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    PICARD_ADDORREPLACEREADGROUPS ( input )
}
