#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAKTA_BAKTA } from '../../../../../modules/nf-core/bakta/bakta/main.nf'

workflow test_bakta_bakta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BAKTA_BAKTA ( input, [], [], [] )
}
