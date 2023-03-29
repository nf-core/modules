#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BACKSUB } from '../../../../modules/nf-core/backsub/main.nf'

workflow test_backsub {

    image = [
        [ id:'test' ], // meta map
        file("/workspace/modules/test_backsub/recyzed_test.ome.tif", checkIfExists: true)
    ]

    markerfile = [
        [ id:'test' ], // meta map
        file("/workspace/modules/test_backsub/markers.csv", checkIfExists: true)
    ]

    BACKSUB ( image, markerfile )
}
