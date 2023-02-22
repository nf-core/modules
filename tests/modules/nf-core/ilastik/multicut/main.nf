#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILASTIK_MULTICUT } from '../../../../../modules/nf-core/ilastik/multicut/main.nf'

workflow test_ilastik_multicut {

    input = [
        [ id:'probabilities' ], // meta map
        file(params.test_data['spatialomics']['h5']['plant_wga'], checkIfExists: true)
    ]
    ilp = [file(params.test_data['spatialomics']['ilp']['plant_wga_multicut'], checkIfExists: true)]
    probs = [file(params.test_data['spatialomics']['h5']['plant_wga_prob'], checkIfExists: true)]

    ILASTIK_MULTICUT ( input, ilp, probs )
}
