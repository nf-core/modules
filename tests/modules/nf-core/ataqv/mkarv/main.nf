#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATAQV_ATAQV } from '../../../../../modules/nf-core/ataqv/ataqv/main.nf'
include { ATAQV_MKARV } from '../../../../../modules/nf-core/ataqv/mkarv/main.nf'

workflow test_ataqv_mkarv {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        [],
        []
    ]

    ATAQV_ATAQV ( input, 'human', '', [], [], [] )
    ATAQV_MKARV ( ATAQV_ATAQV.out.json.collect{ it[1]} )
}
