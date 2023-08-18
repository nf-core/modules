#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STARAMR_SEARCH } from '../../../../../modules/nf-core/staramr/search/main.nf'

workflow test_staramr_search {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    STARAMR_SEARCH ( input )
}

