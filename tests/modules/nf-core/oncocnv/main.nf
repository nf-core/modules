#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONCOCNV } from '../../../../modules/nf-core/oncocnv/main.nf'

workflow test_oncocnv {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    ONCOCNV ( input )
}
