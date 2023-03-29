#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCMICRO_SCIMAP } from '../../../../../modules/nf-core/mcmicro/scimap/main.nf'

workflow test_mcmicro_scimap {

    feature_table = [
        [ id:'test' ], // meta map
        file(params.test_data['imaging']['downstream']['cell_feature_array'], checkIfExists: true)
    ]

    MCMICRO_SCIMAP ( feature_table )
}
