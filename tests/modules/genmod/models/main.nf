#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_MODELS } from '../../../../modules/genmod/models/main.nf'

workflow test_genmod_models {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GENMOD_MODELS ( input )
}
