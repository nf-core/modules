#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SISTR } from '../../../modules/sistr/main.nf'

workflow test_sistr {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true) 
    ]

    SISTR ( input )
}
