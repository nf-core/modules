#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPCELL } from '../../../../modules/nf-core/deepcell/main.nf'

workflow test_deepcell {

    img = [
        [ id:'test_img'], // meta map
        file("/workspace/modules/to_test/Snap-1975.ome.tif", checkIfExists: true)
    ]

    DEEPCELL ( img )
}
