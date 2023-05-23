#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISTOGRAM_HISTOGRAM } from '../../../../../modules/nf-core/histogram/histogram/main.nf'

workflow test_histogram_histogram {
    
    input = [
        [ id:'test' ], // meta map
        "abc123",
        file("/home/ubuntu/abc123_all_stats.txt", checkIfExists: true),
        file("/home/ubuntu/abc123_no_sec_supp_stats.txt", checkIfExists: true)
    ]

    HISTOGRAM_HISTOGRAM ( input )
}
