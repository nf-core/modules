#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAPPA_EXAMINEGRAFT } from '../../../../../modules/nf-core/gappa/examinegraft/main.nf'

workflow test_gappa_examinegraft {
    
    input = [
        [ id:'test' ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/gappa/epa_result.jplace.gz', checkIfExists: true)
    ]

    GAPPA_EXAMINEGRAFT ( input )
}
