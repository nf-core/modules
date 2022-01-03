#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EMMTYPER } from '../../../modules/emmtyper/main.nf'

workflow test_emmtyper {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    EMMTYPER ( input )
}
