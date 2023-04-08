#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINDAGAP_DUPLICATEFINDER } from '../../../../../modules/nf-core/mindagap/duplicatefinder/main.nf'

workflow test_mindagap_duplicatefinder {

    input = [
        [ id:'test' ], // meta map
        file("/workspace/test_data/test_spot_table.tsv", checkIfExists: true)
    ]

    MINDAGAP_DUPLICATEFINDER ( input )
}
