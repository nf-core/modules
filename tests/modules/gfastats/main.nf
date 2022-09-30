#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFASTATS } from '../../../modules/gfastats/main.nf'

workflow test_gfastats {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GFASTATS ( input )
}
