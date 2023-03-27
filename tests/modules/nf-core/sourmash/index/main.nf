#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_INDEX } from '../../../../../modules/nf-core/sourmash/index/main.nf'

workflow test_sourmash_index {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SOURMASH_INDEX ( input )
}
