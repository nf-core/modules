#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PBCCS } from '../../../modules/pbccs/main.nf' addParams( options: [args:'--min-rq 0.9'] )

workflow test_pbccs {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['alz'],    checkIfExists: true),
        file(params.test_data['homo_sapiens']['pacbio']['alzpbi'], checkIfExists: true),
        2,
        3
    ]

    PBCCS ( input )
}
