#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_BEDTOINTERVALLIST } from '../../../../../modules/nf-core/picard/bedtointervallist/main.nf'

workflow test_picard_bedtointervallist {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
    ]
    dict = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    ]

    PICARD_BEDTOINTERVALLIST ( input, dict, [] )
}
