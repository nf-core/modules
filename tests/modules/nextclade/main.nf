#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE } from '../../../modules/nextclade/main.nf'

workflow test_nextclade {
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    
    NEXTCLADE ( input )
}
