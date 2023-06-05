#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { HAMRONIZATION_AMRFINDERPLUS } from '../../../../../modules/nf-core/hamronization/amrfinderplus/main.nf'

workflow test_hamronization_amrfinderplus {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/amrfinderplus/test_output.tsv", checkIfExists: true)
    ]

    HAMRONIZATION_AMRFINDERPLUS ( input, 'tsv', '3.10.30', '2022-05-26.1' )
}
