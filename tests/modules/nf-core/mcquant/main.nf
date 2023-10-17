#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCQUANT } from '../../../../modules/nf-core/mcquant/main.nf'

workflow test_mcquant {

    image = [
        [ id:'image' ], // meta map
        file(params.test_data['imaging']['quantification']['image'], checkIfExists: true)
    ]

    mask = [
        [ id:'mask' ], // meta map
        file(params.test_data['imaging']['quantification']['mask'], checkIfExists: true)
    ]

    markerfile = [
        [ id:'marker' ], // meta map
        file(params.test_data['imaging']['quantification']['markers'], checkIfExists: true)
    ]

    MCQUANT ( image, mask, markerfile )
}
