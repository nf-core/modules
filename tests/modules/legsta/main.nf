#!/usr/bin/env nextflow



include { LEGSTA } from '../../../modules/legsta/main.nf'

workflow test_legsta {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 
    ]

    LEGSTA ( input )
}
