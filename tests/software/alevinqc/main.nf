#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALEVINQC } from '../../../software/alevinqc/main.nf' addParams( options: [:] )

workflow test_alevinqc {

    input = [ [ id:'test'], // meta map
              file("${launchDir}/tests/data/delete_me/alevin_example_v0.14", checkIfExists: true)]

    ALEVINQC ( input )
}
