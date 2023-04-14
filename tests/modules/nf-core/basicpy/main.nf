#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BASICPY } from '../../../../modules/nf-core/basicpy/main.nf'

workflow test_basicpy {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("/workspace/test_data/stack_nocorr.ome.tiff", checkIfExists: true)
    ]

    BASICPY ( input )
}
