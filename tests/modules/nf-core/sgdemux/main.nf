#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SGDEMUX } from '../../../../modules/nf-core/sgdemux/main.nf'
include { UNTAR   } from '../../../../modules/nf-core/untar/main.nf'

workflow test_sgdemux {
    
    input = Channel.value([
        [ id: 'sim-data' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/sim-data/out.sample_meta.csv", checkIfExists: true)
    ])
    ch_input = input.join(
        UNTAR ( [
            [ id:'sim-data' ],
            file("https://github.com/nf-core/test-datasets/blob/demultiplex/testdata/sim-data/fastq.tar.gz?raw=true", checkIfExists: true)
        ]).untar
    )
    SGDEMUX ( ch_input )
}
