#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_PARSE } from '../../../../../modules/nf-core/pairtools/parse/main.nf'

workflow test_pairtools_parse {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['generic']['pairtools']['mock.sam'], checkIfExists: true) ]
    sizes = file(params.test_data['generic']['pairtools']['mock.chrom.sizes'], checkIfExists:true)

    PAIRTOOLS_PARSE ( input, sizes )
}
