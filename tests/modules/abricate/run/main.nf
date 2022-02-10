#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABRICATE_RUN } from '../../../../modules/abricate/run/main.nf'

workflow test_abricate_run {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    ABRICATE_RUN ( input )
}
