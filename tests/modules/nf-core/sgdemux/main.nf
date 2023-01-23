#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SGDEMUX } from '../../../../modules/nf-core/sgdemux/main.nf'

workflow test_sgdemux {
    
    input = Channel.value([
        [ id: 'sim-data' ], // meta map
        file("/Users/swhite/projects/SINGULAR/singular-demux/simulated/out.sample_meta.csv", checkIfExists: true),
        file("/Users/swhite/projects/SINGULAR/singular-demux/simulated/", checkIfExists: true)
    ])

    input.view()
    SGDEMUX ( input )
}
