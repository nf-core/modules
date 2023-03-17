#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCQUANT } from '../../../../modules/nf-core/mcquant/main.nf'

workflow test_mcquant {

    input = [
        [ id:'test' ], // meta map
        file(), checkIfExists: true)
    ]

    MCQUANT ( input )
}
