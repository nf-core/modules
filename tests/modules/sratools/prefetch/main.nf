#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRATOOLS_PREFETCH } from '../../../../modules/sratools/prefetch/main.nf'

workflow test_sratools_prefetch {

    input = [
        [ id:'test', single_end:false ], // meta map
        'ERR2815334'
    ]

    SRATOOLS_PREFETCH ( input )
}
