#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPCELL } from '../../../../modules/nf-core/deepcell/main.nf'

workflow test_deepcell {

    img = [
        [ id:'test_img'], // meta map
        file("/workspace/modules/to_test/cycif_tonsil_small.ome.tif", checkIfExists: true)
    ]
    DEEPCELL ( img )
}
