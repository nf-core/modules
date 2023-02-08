#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINDAGAP } from '../../../../modules/nf-core/mindagap/main.nf'

workflow test_mindagap {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('/Users/florian_wuennemann/1_Projects/nf_core/test-datasets/data/spatialomics/tiff/mindagap.mouse_heart.wga.tiff', checkIfExists: true)
    ]

    MINDAGAP ( input )
}
