#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAMRONIZATION_DEEPARG } from '../../../../modules/hamronization/deeparg/main.nf'

workflow test_hamronization_deeparg {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_mapping_potential_arg'], checkIfExists: true),
    ]

    HAMRONIZATION_DEEPARG ( input, 'tsv', '1.0.2', '2'  )
}
