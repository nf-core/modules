#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TAILFINDR } from '../../../modules/tailfindr/main.nf'

workflow test_tailfindr {
    
    input = [
        [ id: 'test' ], // meta map
        file(params.test_data, checkIfExists: true)
    ]

    TAILFINDR ( input )
}
