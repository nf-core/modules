#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BASICPY } from '../../../../modules/nf-core/basicpy/main.nf'

workflow test_basicpy {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['imaging']['ome-tiff']['cycif_tonsil_cycle1'], checkIfExists: true)
    ]

    BASICPY ( input )
}
