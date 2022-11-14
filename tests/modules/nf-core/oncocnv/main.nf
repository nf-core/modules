#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONCOCNV } from '../../../../modules/nf-core/oncocnv/main.nf'

workflow test_oncocnv {
    
    normal = [
        "normal",
        [
            file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
        ],
        [
            file(params.test_data['sarscov2']['illumina']['test_paired_end_bai'], checkIfExists: true)
        ],
    ]

    ONCOCNV ( normal )
}
