#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_DUMP } from '../../../../software/cooler/dump/main.nf' addParams( options: [:] )

workflow test_cooler_dump {

    input = [ [ id:'test', bin:16 ], // meta map
              file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/toy.asymm.16.cool", checkIfExists: true) ]

    COOLER_DUMP ( input )
}
