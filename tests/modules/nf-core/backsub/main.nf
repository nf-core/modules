#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BACKSUB } from '../../../../modules/nf-core/backsub/main.nf'

workflow test_backsub {

    image = [
        [ id:'test' ], // meta map
        file(params.test_data['imaging']['background_subtraction']['image'], checkIfExists: true)
    ]

    markerfile = [
        [ id:'test' ], // meta map
        file(params.test_data['imaging']['background_subtraction']['marker'], checkIfExists: true)
    ]

    BACKSUB ( image, markerfile )
}
