#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_MERGE } from '../../../../modules/cooler/merge/main.nf' addParams( options: [:] )

workflow test_cooler_merge {

    input = [ [ id:'test' ], // meta map
              [ file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/toy.symm.upper.2.cool", checkIfExists: true) ] ]

    COOLER_MERGE ( input )
}
