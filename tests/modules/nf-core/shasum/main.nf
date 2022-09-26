#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHASUM } from '../../../modules/shasum/main.nf'

workflow test_shasum {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SHASUM ( input )
}
