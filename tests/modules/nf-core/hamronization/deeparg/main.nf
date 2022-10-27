#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { HAMRONIZATION_DEEPARG } from "$moduleDir/modules/nf-core/hamronization/deeparg/main.nf"

workflow test_hamronization_deeparg {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['hamronization']['genome_mapping_potential_arg'], checkIfExists: true),
    ]

    HAMRONIZATION_DEEPARG ( input, 'tsv', '1.0.2', '2'  )
}
