#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_MERGE } from '../../../../software/cooler/merge/main.nf' addParams( options: [:] )

workflow test_cooler_merge {

    input = Channel.of(file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/toy.symm.upper.2.cool", checkIfExists: true))
                   .collect()
                   .map{ [[ id:'test', bin:2 ], it] }

    COOLER_MERGE ( input )
}
