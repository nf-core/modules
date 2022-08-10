#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TEST } from '../../../modules/test/main.nf'

workflow test_test {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    TEST ( input )
}
