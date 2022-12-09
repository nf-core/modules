#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_PURGEDUPS } from '../../../../../modules/nf-core/purgedups/purgedups/main.nf'

workflow test_purgedups_purgedups {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PURGEDUPS_PURGEDUPS ( input )
}
