#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHANTOMPEAKQUALTOOLS } from '../../../modules/phantompeakqualtools/main.nf'

workflow test_phantompeakqualtools_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
    ]

    PHANTOMPEAKQUALTOOLS ( input )
}

workflow test_phantompeakqualtools_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PHANTOMPEAKQUALTOOLS ( input )
}
