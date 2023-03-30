#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SCIMAP_MCMICRO } from '../../../../../modules/nf-core/scimap/mcmicro/main.nf'

workflow test_scimap_mcmicro {

    feature_table = [
        [ id:'test' ], // meta map
        file(params.test_data['imaging']['downstream']['cell_feature_array'], checkIfExists: true)
    ]

    SCIMAP_MCMICRO ( feature_table )
}
