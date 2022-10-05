#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIBAMSUMMARY } from '../../../modules/multibamsummary/main.nf'

workflow test_multibamsummary {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
        // add test_paired_end_sorted_bam and test_paired_end_sorted_bam_bai?
    ]

    MULTIBAMSUMMARY ( input )
}
