#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MIRANDA } from '../../../../modules/nf-core/miranda/main.nf'

workflow test_miranda {

    input = [
        [ id:'cel_N2_1' ],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/miranda/cel_N2_1.fa")
    ]

    mirbase = [
        file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/miranda/mature.fa")
    ]

    MIRANDA ( input, mirbase )
}
