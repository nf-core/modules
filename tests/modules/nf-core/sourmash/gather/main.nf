#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_SKETCH } from '../../../../../modules/nf-core/sourmash/sketch/main.nf'
include { SOURMASH_GATHER } from '../../../../../modules/nf-core/sourmash/gather/main.nf'

workflow test_sourmash_gather {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    //// Database with reference signatures
    //database = [
    //    [ id:'test', single_end:false ], // meta map
    //    file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    //]

    // For the sake of the test simplicity,
    // the same signatures will be used as an input and as database
    // However, there will be no unassigned signatures in the results

    SOURMASH_SKETCH ( input )
    SOURMASH_GATHER ( SOURMASH_SKETCH.out.signatures, SOURMASH_SKETCH.out.signatures )
}
