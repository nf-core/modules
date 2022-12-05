#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_PBCSTAT } from '../../../../../modules/nf-core/purgedups/pbcstat/main.nf'

workflow test_purgedups_pbcstat {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PURGEDUPS_PBCSTAT ( input )
}
