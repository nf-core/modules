#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASSEMBLYSCAN } from '../../../modules/assemblyscan/main.nf'

workflow test_assemblyscan {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    ASSEMBLYSCAN ( input )
}
