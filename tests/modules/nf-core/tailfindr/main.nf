#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TAILFINDR } from '../../../../modules/nf-core/tailfindr/main.nf'

workflow test_tailfindr {
    
    input = [
        [ id: 'test' ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/tailfindr/test.fast5', checkIfExists: true)
    ]

    TAILFINDR ( input )
}
