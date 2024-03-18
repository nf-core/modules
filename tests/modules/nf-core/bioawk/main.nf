#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BIOAWK } from '../../../../modules/nf-core/bioawk/main.nf'

workflow test_bioawk {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BIOAWK ( input )
}
