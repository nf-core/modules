#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { COOLER_DUMP } from "$moduleDir/modules/nf-core/cooler/dump/main.nf"

workflow test_cooler_dump {

    input = [ [ id:'test' ], // meta map
              file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/toy.asymm.16.cool", checkIfExists: true) ]

    COOLER_DUMP ( input, [:] )
}
