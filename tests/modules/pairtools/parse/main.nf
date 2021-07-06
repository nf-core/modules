#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_PARSE } from '../../../../software/pairtools/parse/main.nf' addParams( options: ['suffix':'.raw'] )

workflow test_pairtools_parse {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.sam", checkIfExists: true) ]
    sizes = file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.chrom.sizes", checkIfExists:true)

    PAIRTOOLS_PARSE ( input, sizes )
}
