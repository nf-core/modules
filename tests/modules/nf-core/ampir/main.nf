#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPIR } from '../../../../modules/nf-core/ampir/main.nf'

workflow test_ampir {

    fasta = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['proteome_fasta'], checkIfExists: true),
    ]

    model = "precursor"

    min_length = 10

    min_probability = "0.7"

    AMPIR ( fasta, model, min_length, min_probability )
}
