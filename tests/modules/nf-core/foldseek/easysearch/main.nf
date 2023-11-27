#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOLDSEEK_CREATEDB   } from '../../../../../modules/nf-core/foldseek/createdb/main.nf'
include { FOLDSEEK_EASYSEARCH } from '../../../../../modules/nf-core/foldseek/easysearch/main.nf'

workflow test_foldseek_easysearch {

    input = [
        [ id:'test_db' ],
        [ file(params.test_data['proteomics']['pdb']['tim1_pdb'], checkIfExists: true) ]
    ]

    input2 = [
        [ id:'test_search' ],
        [ file(params.test_data['proteomics']['pdb']['tim8_pdb'], checkIfExists: true) ]
    ]

    FOLDSEEK_CREATEDB ( input )
    FOLDSEEK_EASYSEARCH ( input2, FOLDSEEK_CREATEDB.out.db )
}
