#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPATYPER } from '../../../modules/spatyper/main.nf'
include { SPATYPER as SPATYPER_ENRICH } from '../../../modules/spatyper/main.nf'

workflow test_spatyper {
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    repeats = []
    repeat_order = []

    SPATYPER ( input, repeats, repeat_order )
}

workflow test_spatyper_enrich {
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    repeats = []
    repeat_order = []

    SPATYPER_ENRICH ( input, repeats, repeat_order )
}
