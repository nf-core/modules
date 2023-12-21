#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GPROFILER2_GOST } from '../../../../../modules/nf-core/gprofiler2/gost/main.nf'

workflow test_gprofiler2_gost {
    contrasts = [ [ id:'test', reference:'r', target:'t' ], 'test', 'r', 't' ]
    input = [
        [ id:'test_r_t', reference:'r', target:'t' ], // meta map
        file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true)
    ]

    GPROFILER2_GOST (
        input,
        [],
        []
    )
}

workflow test_gprofiler2_gost_backgroundmatrix {
    contrasts = [ [ id:'test', reference:'r', target:'t' ], 'test', 'r', 't' ]
    input = [
        [ id:'test_r_t', reference:'r', target:'t' ], // meta map
        file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true)
    ]
    ch_background = Channel.from(file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true))

    GPROFILER2_GOST (
        input,
        [],
        ch_background
    )
}
