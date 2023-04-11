#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLPOSE } from '../../../../modules/nf-core/cellpose/main.nf'

workflow test_cellpose {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['imaging']['segmentation']['image'], checkIfExists: true)
    ]

    CELLPOSE ( input, [] )
}
