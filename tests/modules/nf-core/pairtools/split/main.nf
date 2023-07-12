#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_SPLIT } from '../../../../../modules/nf-core/pairtools/split/main.nf'

workflow test_pairtools_dedup {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['generic']['pairtools']['mock_4dedup_pairsam'], checkIfExists: true) ]

    PAIRTOOLS_SPLIT ( input )
}

