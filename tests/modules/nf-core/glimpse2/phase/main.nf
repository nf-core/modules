#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE2_PHASE } from '../../../../../modules/nf-core/glimpse2/phase/main.nf'

workflow test_glimpse2_phase {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GLIMPSE2_PHASE ( input )
}
