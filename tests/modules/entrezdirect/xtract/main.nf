#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_XTRACT } from '../../../../modules/entrezdirect/xtract/main.nf'

workflow test_entrezdirect_xtract {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    ENTREZDIRECT_XTRACT ( input )
}
