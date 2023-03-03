#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_STATS } from '../../../../../modules/nf-core/pairtools/stats/main.nf'

workflow test_pairtools_stats {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.pairsam', checkIfExists: true)
    ]

    PAIRTOOLS_STATS ( input )
}
