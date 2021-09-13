#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPATYPER } from '../../../modules/spatyper/main.nf' addParams( options: [:] )
include { SPATYPER as SPATYPER_ENRICH } from '../../../modules/spatyper/main.nf' addParams( options: [args: '--do_enrich'] )

workflow test_spatyper {
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    SPATYPER ( input )
}

workflow test_spatyper_enrich {
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    SPATYPER_ENRICH ( input )
}
