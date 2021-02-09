#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../software/untar/main.nf' addParams( options: [:] )

workflow test_untar {

    def input = []
    input = [ file("${launchDir}/tests/data/tar_gz/test.txt.tar.gz", checkIfExists: true) ]

    UNTAR ( input )
}
