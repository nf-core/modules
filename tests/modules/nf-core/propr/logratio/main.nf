#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROPR_LOGRATIO } from '../../../../../modules/nf-core/propr/logratio/main.nf'

workflow test_propr_logratio {
    
    input = [
        [ id:'test.clr.NA.NA', gene_id_col:'gene_id' ], // meta map
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true),
        'clr',
        null,
        null
    ]

    PROPR_LOGRATIO ( input )
}
