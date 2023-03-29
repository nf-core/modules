#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLPOSE } from '../../../../modules/nf-core/cellpose/main.nf'

workflow test_cellpose {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("/workspace/modules/modules/nf-core/cellpose/cycif_tonsil_small.ome.tif", checkIfExists: true)
    ]

    CELLPOSE ( input )
}
