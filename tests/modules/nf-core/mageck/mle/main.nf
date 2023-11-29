#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAGECK_MLE } from '../../../../../modules/nf-core/mageck/mle/main.nf'

workflow test_mageck_mle {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['mus_musculus']['csv']['count_table'], checkIfExists: true)
    ]
    design_matrix = file(params.test_data['mus_musculus']['txt']['design_matrix'], checkIfExists: true)
    MAGECK_MLE ( input, design_matrix )
}
