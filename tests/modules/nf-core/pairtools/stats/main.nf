#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_STATS } from '../../../../../modules/nf-core/pairtools/stats/main.nf'

workflow test_pairtools_stats {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['generic']['cooler']['test_pairs_pair'], checkIfExists: true)
    ]

    PAIRTOOLS_STATS ( input )
}
