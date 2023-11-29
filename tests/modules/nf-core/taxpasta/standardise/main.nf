#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TAXPASTA_STANDARDISE } from '../../../../../modules/nf-core/taxpasta/standardise/main.nf'

workflow test_taxpasta_standardise {

    input = [
        [ id:'test', single_end:false, tool:'kraken2' ], // meta map
        file(params.test_data['sarscov2']['metagenome']['kraken_report'], checkIfExists: true)
    ]

    TAXPASTA_STANDARDISE(input, [])
}
