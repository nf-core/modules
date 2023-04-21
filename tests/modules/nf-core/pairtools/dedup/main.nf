#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_DEDUP } from '../../../../../modules/nf-core/pairtools/dedup/main.nf'

workflow test_pairtools_dedup {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['generic']['pairtools']['mock.4dedup.pairsam'], checkIfExists: true) ]

    PAIRTOOLS_DEDUP ( input )
}
