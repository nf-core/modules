#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_GETSEQS } from '../../../../../modules/nf-core/purgedups/getseqs/main.nf'

workflow test_purgedups_getseqs {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PURGEDUPS_GETSEQS ( input )
}
