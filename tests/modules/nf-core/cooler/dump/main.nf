#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_DUMP } from '../../../../../modules/nf-core/cooler/dump/main.nf'

workflow test_cooler_dump {

    input = [ [ id:'test' ], // meta map
              file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/toy.asymm.16.cool", checkIfExists: true),
              [:] // resolution if any
    ]

    COOLER_DUMP ( input )
}
