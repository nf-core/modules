#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_ESLALIMASK } from '../../../../modules/hmmer/eslalimask/main.nf'

workflow test_hmmer_eslalimask {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    HMMER_ESLALIMASK ( input )
}
