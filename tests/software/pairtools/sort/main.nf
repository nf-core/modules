#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_SORT } from '../../../../software/pairtools/sort/main.nf' addParams( options: ['suffix':'.sorted'] )

workflow test_pairtools_sort {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.pairsam", checkIfExists: true) ]

    PAIRTOOLS_SORT ( input )
}
