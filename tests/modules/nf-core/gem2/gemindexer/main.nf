#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GEM2_GEMINDEXER } from '../../../../../modules/nf-core/gem2/gemindexer/main.nf'

workflow test_gem2_gemindexer {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    GEM2_GEMINDEXER ( input )
}
