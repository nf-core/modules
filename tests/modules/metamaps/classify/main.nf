#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METAMAPS_CLASSIFY } from '../../../../modules/metamaps/classify/main.nf'

workflow test_metamaps_classify {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    METAMAPS_CLASSIFY ( input )
}
