#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUBBINS } from '../../../software/gubbins/main.nf' addParams( options: [:] )

workflow test_gubbins { 
    input = file("${launchDir}/tests/data/alignment/test.fas", checkIfExists: true)

    GUBBINS ( input )
}