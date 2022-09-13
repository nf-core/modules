#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BASES2FASTQ } from '../../../modules/bases2fastq/main.nf'

workflow test_bases2fastq {
    
    input = [
        [ id:'sim-data' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/sim-data/RunManifest.csv", checkIfExists: true),
        file("s3://element-public-data/sim-data/", checkIfExists: true)
    ]

    BASES2FASTQ ( input )
}
