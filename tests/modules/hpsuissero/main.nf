#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HPSUISSERO } from '../../../modules/hpsuissero/main.nf'

workflow test_hpsuissero {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true) 
    ]

    HPSUISSERO ( input )
}
