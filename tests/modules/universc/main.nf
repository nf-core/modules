#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNIVERSC } from '../../../modules/universc/main.nf'

workflow test_universc {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    UNIVERSC ( input )
}
