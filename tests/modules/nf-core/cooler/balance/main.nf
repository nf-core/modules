#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_BALANCE } from '../../../../../modules/nf-core/cooler/balance/main.nf'
include { COOLER_DUMP    } from '../../../../../modules/nf-core/cooler/dump/main.nf'

workflow test_cooler_balance {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['generic']['cooler']['test_merge_cool'], checkIfExists: true),
        ''
    ]

    COOLER_BALANCE ( input )
    COOLER_DUMP ( COOLER_BALANCE.out.cool.combine([[:]]) )
}
