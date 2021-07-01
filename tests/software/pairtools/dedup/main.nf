#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_DEDUP } from '../../../../software/pairtools/dedup/main.nf' addParams( options: ['suffix':'.dedup'] )

workflow test_pairtools_dedup {

    input = [ [ id:'test', single_end:false ], // meta map
                file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.4dedup.pairsam", checkIfExists: true) ]

    PAIRTOOLS_DEDUP ( input )
}
