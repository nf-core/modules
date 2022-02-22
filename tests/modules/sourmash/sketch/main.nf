#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_SKETCH } from '../../../../modules/sourmash/sketch/main.nf'

workflow test_sourmash_sketch {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    SOURMASH_SKETCH ( input )
}
