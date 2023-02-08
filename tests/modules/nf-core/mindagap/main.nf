#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINDAGAP } from '../../../../modules/nf-core/mindagap/main.nf'

workflow test_mindagap {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['spatialomics']['tiff']['mouse_heart_wga_tiff'], checkIfExists: true)
    ]

    MINDAGAP ( input )
}
