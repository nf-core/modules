#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABRICATE_ABRICATE } from '../../../../modules/abricate/abricate/main.nf'

workflow test_abricate_abricate {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    ABRICATE_ABRICATE ( input )
}
