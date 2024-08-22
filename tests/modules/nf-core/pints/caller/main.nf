#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PINTS_CALLER } from '../../../../../modules/nf-core/pints/caller/main.nf'

// This tests a single bam input, and
// ttps://github.com/hyulab/PINTS/issues/12
workflow test_pints_caller_empty_results {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    PINTS_CALLER ( input )
}

// TODO Test single bigwig input
// TODO Test multiple bigwig input
