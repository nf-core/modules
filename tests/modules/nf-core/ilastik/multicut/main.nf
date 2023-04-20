#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILASTIK_MULTICUT } from '../../../../../modules/nf-core/ilastik/multicut/main.nf'

workflow test_ilastik_multicut {

    input = [
        [ id:'image' ], // meta map
        file(params.test_data['imaging']['h5']['plant_wga'], checkIfExists: true)
    ]
    ilp =   [
        [id:'ilastik_model'],
        file(params.test_data['imaging']['ilp']['plant_wga_multicut'], checkIfExists: true)
    ]
    probs = [
        [id:'probability_maps'],
        file(params.test_data['imaging']['h5']['plant_wga_prob'], checkIfExists: true)
    ]

    ILASTIK_MULTICUT ( input, ilp, probs )
}
