#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_SCORE } from '../../../../modules/genmod/score/main.nf'

workflow test_genmod_score {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GENMOD_SCORE ( input )
}
