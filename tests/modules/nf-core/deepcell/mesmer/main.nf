#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPCELL_MESMER } from '../../../../../modules/nf-core/deepcell/mesmer/main.nf'

workflow test_deepcell_mesmer {

    img = [
        [ id:'test_img'], // meta map
        file(params.test_data['imaging']['segmentation']['image'], checkIfExists: true)
    ]

    DEEPCELL_MESMER ( img,  [[:],[]])
}
