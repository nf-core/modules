#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTINSERTSIZEMETRICS } from '../../../../../modules/nf-core/picard/collectinsertsizemetrics/main.nf'

workflow test_picard_collectinsertsizemetrics {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PICARD_COLLECTINSERTSIZEMETRICS ( input )
}
