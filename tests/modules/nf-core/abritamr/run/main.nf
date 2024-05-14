#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABRITAMR_RUN } from '../../../../../modules/nf-core/abritamr/run/main.nf'

workflow test_abritamr_run {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    ABRITAMR_RUN ( input )
}
