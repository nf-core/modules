#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CRISPRCLEANR_NORMALIZE } from '../../../../../modules/nf-core/crisprcleanr/normalize/main.nf'
workflow test_crisprcleanr_normalize {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['crispr_functional_counts'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['crispr_functional_library'], checkIfExists: true)
    ]

    CRISPRCLEANR_NORMALIZE ( input, 30, 5)
}
