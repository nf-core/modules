#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_SKETCH as SOURMASH_SKETCH1 } from '../../../../../modules/nf-core/sourmash/sketch/main.nf'
include { SOURMASH_SKETCH as SOURMASH_SKETCH2 } from '../../../../../modules/nf-core/sourmash/sketch/main.nf'
include { SOURMASH_GATHER } from '../../../../../modules/nf-core/sourmash/gather/main.nf'

workflow test_sourmash_gather {

    // For the sake of the test simplicity,
    // the same signatures will be used as an input and as database

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    database = [
        [ id:'db', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    SOURMASH_SKETCH1 ( input )
    SOURMASH_SKETCH2 ( database )

    SOURMASH_GATHER ( SOURMASH_SKETCH1.out.signatures, SOURMASH_SKETCH2.out.signatures.map{it[1]} )
}
