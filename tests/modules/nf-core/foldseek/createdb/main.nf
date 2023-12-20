#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDSEEK_CREATEDB } from '../../../../../modules/nf-core/foldseek/createdb/main.nf'

workflow test_foldseek_createdb {

    input = [
        [ id:'test' ],
        file(params.test_data['proteomics']['pdb']['tim1_pdb'], checkIfExists: true)
    ]

    FOLDSEEK_CREATEDB ( input )
}
