#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROKKA } from '../../../modules/prokka/main.nf'

workflow test_prokka {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    
    PROKKA ( input, [], [] )
}
