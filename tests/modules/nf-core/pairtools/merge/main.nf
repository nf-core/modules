#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_MERGE } from '../../../../../modules/nf-core/pairtools/merge/main.nf'

workflow test_pairtools_merge {
    
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['generic']['pairtools']['mock_4dedup_pairsam'], checkIfExists: true) ]

    PAIRTOOLS_MERGE ( input )
}
