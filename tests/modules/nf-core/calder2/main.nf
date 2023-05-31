#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CALDER2 } from '../../../../modules/nf-core/calder2/main.nf'


def TEST_COOL_PATH = "https://raw.githubusercontent.com/CSOgroup/CALDER2/main/tests/testthat/data/test.cool"
def TEST_MCOOL_PATH = "https://raw.githubusercontent.com/CSOgroup/CALDER2/main/tests/testthat/data/test.mcool"


workflow test_calder2_cool {
    
    input = [ [ id:'test' ], //meta map
              file(TEST_COOL_PATH, checkIfExists: true) ]

    CALDER2 ( input, [:] )
}

workflow test_calder2_mcool {
    input = [ [ id:'test' ], //meta map
              file(TEST_MCOOL_PATH, checkIfExists: true) ]

    CALDER2 ( input, 100000 )
}