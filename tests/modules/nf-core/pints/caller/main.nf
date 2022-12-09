#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PINTS_CALLER } from '../../../../../modules/nf-core/pints/caller/main.nf'

workflow test_pints_caller {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        // FIXME Fails if it doesn't find any signals
        [file("https://raw.githubusercontent.com/Kraus-Lab/groHMM/master/inst/extdata/S0mR1.bam", checkIfExists: true),
         file("https://raw.githubusercontent.com/Kraus-Lab/groHMM/master/inst/extdata/S40mR1.bam", checkIfExists: true)]
    ]

    PINTS_CALLER ( input )
}

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
