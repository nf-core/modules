#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_SELECT } from '../../../../modules/pairtools/select/main.nf'

workflow test_pairtools_select {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.pairsam", checkIfExists: true) ]

    PAIRTOOLS_SELECT ( input )
}
