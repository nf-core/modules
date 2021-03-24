#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../software/untar/main.nf' addParams( options: [:] )

workflow test_untar {
    input = file("${launchDir}/tests/data/generic/test.txt.tar.gz", checkIfExists: true)

    UNTAR ( input )
}
