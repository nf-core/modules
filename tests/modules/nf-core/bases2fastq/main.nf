#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BASES2FASTQ } from '../../../../modules/nf-core/bases2fastq/main.nf'
include { UNTAR      } from '../../../../modules/nf-core/untar/main.nf'

workflow test_bases2fastq {
    
    input = Channel.value([
        [ id:'sim-data' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/sim-data/RunManifest.csv", checkIfExists: true),
    ])

    ch_input = input.join(
        UNTAR ( [
            [ id:'sim-data' ],
            file("https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/sim-data/sim-data.tar.gz", checkIfExists: true)
        ]).untar
    )

    BASES2FASTQ ( ch_input )
}
