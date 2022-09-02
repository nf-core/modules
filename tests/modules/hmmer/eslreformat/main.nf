#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_ESLREFORMAT } from '../../../../modules/hmmer/eslreformat/main.nf'

workflow test_hmmer_eslreformat {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    HMMER_ESLREFORMAT ( input )
}
