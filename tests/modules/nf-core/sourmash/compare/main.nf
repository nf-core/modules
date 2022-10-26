#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_SKETCH } from '../../../../../modules/nf-core/sourmash/sketch/main.nf'
include { SOURMASH_COMPARE } from '../../../../../modules/nf-core/sourmash/compare/main.nf'

workflow test_sourmash_compare {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    save_numpy_matrix = true
    save_csv = true
    
    SOURMASH_SKETCH ( input )

    SOURMASH_COMPARE (
        SOURMASH_SKETCH.out.signatures.collect { it[1] },
        save_numpy_matrix,
        save_csv
     )
}
