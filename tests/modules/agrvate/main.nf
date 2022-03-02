#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGRVATE } from '../../../modules/agrvate/main.nf'

workflow test_agrvate {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    AGRVATE ( input )
}
