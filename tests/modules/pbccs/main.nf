#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PBCCS } from '../../../modules/pbccs/main.nf' addParams( options: [args:'--min-rq 0.9'] )

workflow test_pbccs {

    input = [
        [ id:'alz' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['alz'],    checkIfExists: true),
        file(params.test_data['homo_sapiens']['pacbio']['alzpbi'], checkIfExists: true)
    ]

    chunk_num = 2
    chunk_on  = 3

    PBCCS ( input, chunk_num, chunk_on )
}
