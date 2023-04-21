#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_SORT } from '../../../../../modules/nf-core/pairtools/sort/main.nf'

workflow test_pairtools_sort {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['generic']['pairtools']['mock.pairsam'], checkIfExists: true) ]

    PAIRTOOLS_SORT ( input )
}
