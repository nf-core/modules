#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SUBTRACT } from '../../../../modules/bedtools/subtract/main.nf'

workflow test_bedtools_subtract {
    input = [
        [ id:'test_subtract' ],
        file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    BEDTOOLS_SUBTRACT ( input )
}
