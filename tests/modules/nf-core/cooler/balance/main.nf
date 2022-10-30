#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_BALANCE } from '../../../../../modules/nf-core/cooler/balance/main.nf'

workflow test_cooler_balance {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['generic']['cooler']['test_merge_cool'], checkIfExists: true),
	''
    ]

    COOLER_BALANCE ( input )
}
