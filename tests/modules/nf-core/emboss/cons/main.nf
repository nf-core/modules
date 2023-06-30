#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EMBOSS_CONS } from '../../../../../modules/nf-core/emboss/cons/main.nf'

workflow test_emboss_cons {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['all_sites_fas'], checkIfExists: true)
    ]

    EMBOSS_CONS ( input )
}
