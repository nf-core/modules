#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../../modules/picard/createsequencedictionary/main.nf'

workflow test_picard_createsequencedictionary {
    
    input = [ 
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 
    ]

    PICARD_CREATESEQUENCEDICTIONARY ( input )
}
