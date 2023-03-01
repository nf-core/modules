#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHAPEIT5_SWITCH } from '../../../../../modules/nf-core/shapeit5/switch/main.nf'

workflow test_shapeit5_switch {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SHAPEIT5_SWITCH ( input )
}
