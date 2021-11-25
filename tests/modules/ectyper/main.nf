#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ECTYPER } from '../../../modules/ectyper/main.nf'

workflow test_ectyper {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ECTYPER ( input )
}
