#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_SKETCH } from '../../../../../modules/nf-core/sourmash/sketch/main.nf'
include { SOURMASH_COMPARE } from '../../../../../modules/nf-core/sourmash/compare/main.nf'

workflow test_sourmash_compare {
    
    input = Channel.fromList([
       [ [ id:'test1', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ],
       [ [ id: 'test2', single_end:false ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    ])

    ch_sketch_for_compare = SOURMASH_SKETCH ( input ).signatures
        .collect { it[1] }
        .map {
            signatures ->
                def meta = [:]
                meta.id = "group1"
                [ meta, signatures ]
        }

    save_numpy_matrix = true
    save_csv = true

    SOURMASH_COMPARE (
        ch_sketch_for_compare,
        [],
        save_numpy_matrix,
        save_csv
     )
    
    SOURMASH_COMPARE.out.matrix
}
