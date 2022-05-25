#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPIR } from '../../../modules/ampir/main.nf'

workflow test_ampir {

    fasta = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['proteome_fasta'], checkIfExists: true),
    ]

    cut_off = "80"

    model = "precursor"

    output_name = "prediction.fasta"

    AMPIR ( fasta, cut_off, model, output_name )
}
