#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { SRATOOLS_PREFETCH } from "$moduleDir/modules/nf-core/sratools/prefetch/main.nf"

workflow test_sratools_prefetch {

    input = [
        [ id:'test', single_end:false ], // meta map
        'DRR000774'
    ]

    SRATOOLS_PREFETCH(input, file(params.test_data['generic']['config']['ncbi_user_settings'], checkIfExists: true))
}
