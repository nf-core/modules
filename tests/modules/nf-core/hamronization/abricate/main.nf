#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAMRONIZATION_ABRICATE } from '../../../../../modules/nf-core/hamronization/abricate/main.nf'

workflow test_hamronization_abricate {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['hamronization']['genome_abricate_tsv'], checkIfExists: true),
    ]

    HAMRONIZATION_ABRICATE ( input, 'tsv', '1.0.1', '2021-Mar-27' )
}
