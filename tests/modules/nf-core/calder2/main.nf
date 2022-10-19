#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CALDER2 } from '../../../../modules/nf-core/calder2/main.nf'


def TEST_COOL_PATH = "https://raw.githubusercontent.com/CSOgroup/CALDER2/main/tests/testthat/data/test.cool"


workflow test_calder2 {
    
    input = [ [ id:'test' ], //meta map
              file(TEST_COOL_PATH, checkIfExists: true) ]

    CALDER2 ( input, [:] )
}
