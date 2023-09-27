#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROPR_LOGRATIO } from '../../../../../modules/nf-core/propr/logratio/main.nf'

workflow test_propr_logratio {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ]

    PROPR_LOGRATIO ( input )
}
