#!/usr/bin/env nextflow



include { SSUISSERO } from '../../../modules/ssuissero/main.nf'

workflow test_ssuissero {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true) 
    ]

    SSUISSERO ( input )
}
