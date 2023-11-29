#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAPPA_EXAMINEASSIGN } from '../../../../../modules/nf-core/gappa/examineassign/main.nf'

workflow test_gappa_examineassign {
    
    input = [
        [ id:'test' ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/gappa/epa_result.jplace.gz', checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/gappa/gappa_taxonomy.tsv', checkIfExists: true)
    ]

    GAPPA_EXAMINEASSIGN ( input )
}
