#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONIFER } from '../../../../modules/nf-core/conifer/main.nf'

workflow test_conifer {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    CONIFER ( input, [] )
}