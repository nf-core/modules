#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_FLIP } from '../../../../modules/pairtools/flip/main.nf'

workflow test_pairtools_flip {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.4flip.pairs", checkIfExists: true) ]
    sizes = file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.chrom.sizes", checkIfExists:true)

    PAIRTOOLS_FLIP ( input, sizes )
}
