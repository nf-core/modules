#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRATOOLS_PREFETCH } from '../../../../../modules/nf-core/sratools/prefetch/main.nf'

workflow test_sratools_prefetch {

    input = [
        [ id:'test', single_end:false ], // meta map
        'DRR000774'
    ]

    SRATOOLS_PREFETCH(input, file(params.test_data['generic']['config']['ncbi_user_settings'], checkIfExists: true))
}

workflow test_sratools_prefetch_with_sralite {

    input = [
        [ id:'test', single_end:false ], // meta map
        'SRR1170046'
    ]

    SRATOOLS_PREFETCH(input, file(params.test_data['generic']['config']['ncbi_user_settings'], checkIfExists: true))
}
