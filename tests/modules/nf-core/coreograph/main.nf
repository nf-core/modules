#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COREOGRAPH } from '../../../../modules/nf-core/coreograph/main.nf'

workflow test_coreograph {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['imaging']['core_detection']['image'], checkIfExists: true)
    ]

    COREOGRAPH ( input )
}
