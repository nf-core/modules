#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCQUANT } from '../../../../modules/nf-core/mcquant/main.nf'

workflow test_mcquant {

    image = [
        [ id:'test' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/mcquant/image.tiff", checkIfExists: true) // TODO : replace with public test data
    ]

    mask = [
        [ id:'test' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/mcquant/mask.tif", checkIfExists: true)  // TODO : replace with public test data
    ]

    markerfile = [
        [ id:'test' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/mcquant/channel_names.csv", checkIfExists: true)  // TODO : replace with public test data
    ]

    MCQUANT ( image, mask, markerfile )
}
