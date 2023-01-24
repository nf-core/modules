#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SGDEMUX } from '../../../../modules/nf-core/sgdemux/main.nf'
include { UNTAR      } from '../../../../modules/nf-core/untar/main.nf'

workflow test_sgdemux {
    
    input = Channel.value([
        [ id: 'sim-data' ], // meta map
        file("/Users/swhite/projects/SINGULAR/singular-demux/simulated/out.sample_meta.csv", checkIfExists: true)
        //file("/Users/swhite/projects/SINGULAR/singular-demux/simulated/", checkIfExists: true)
    ])
    ch_input = input.join(
        UNTAR ( [
            [ id:'sim-data' ],
            // tar ball created with tar -czf fastq.tar.gz *.fastq.gz
            file("/Users/swhite/projects/SINGULAR/singular-demux/simulated/fastq.tar.gz", checkIfExists: true)
        ]).untar
    )
    SGDEMUX ( ch_input )
}
