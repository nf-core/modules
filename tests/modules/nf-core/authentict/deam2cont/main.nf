#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AUTHENTICT_DEAM2CONT } from '../../../../../modules/nf-core/authentict/deam2cont/main.nf'

workflow test_authentict_deam2cont {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    AUTHENTICT_DEAM2CONT ( input, [[], []], [[], []] )
}
