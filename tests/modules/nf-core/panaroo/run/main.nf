#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PANAROO_RUN } from '../../../../modules/panaroo/run/main.nf'

workflow test_panaroo_run {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['test1_gff'], checkIfExists: true),
            file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['test2_gff'], checkIfExists: true),
            file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['test3_gff'], checkIfExists: true)
        ]
    ]

    PANAROO_RUN ( input )
}
