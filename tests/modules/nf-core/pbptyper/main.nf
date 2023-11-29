#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PBPTYPER } from '../../../../modules/nf-core/pbptyper/main.nf'

workflow test_pbptyper {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    PBPTYPER ( input, [] )
}
