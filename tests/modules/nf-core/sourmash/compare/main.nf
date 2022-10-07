#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_COMPARE } from '../../../../../modules/nf-core/sourmash/compare/main.nf'

workflow test_sourmash_compare {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SOURMASH_COMPARE ( input )
}
