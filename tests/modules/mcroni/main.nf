#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCRONI } from '../../../modules/mcroni/main.nf'

workflow test_mcroni {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    MCRONI ( input )
}
