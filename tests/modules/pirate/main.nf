#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PIRATE } from '../../../modules/pirate/main.nf'

workflow test_pirate {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['test1_gff'], checkIfExists: true),
            file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['test2_gff'], checkIfExists: true),
            file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['test3_gff'], checkIfExists: true)
        ]
    ]
    //           [ file("https://github.com/bactopia/bactopia-tests/raw/main/data/reference/gff/GCF_000292685.gff", checkIfExists: true),
    //             file("https://github.com/bactopia/bactopia-tests/raw/main/data/reference/gff/GCF_000298385.gff", checkIfExists: true),
    //             file("https://github.com/bactopia/bactopia-tests/raw/main/data/reference/gff/GCF_002849995.gff", checkIfExists: true) ]
    // ]

    PIRATE ( input )
}
