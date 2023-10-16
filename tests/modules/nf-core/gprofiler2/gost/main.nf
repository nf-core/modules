#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GPROFILER2_GOST } from '../../../../../modules/nf-core/gprofiler2/gost/main.nf'

workflow test_gprofiler2_gost {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true)
    ]

    GPROFILER2_GOST (
        input,
        []
    )
}


workflow test_gprofiler2_gost_backgroundmatrix {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true)
    ]
   /* background = [
        [ id:'test'],
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ] */
    background = Channel.fromPath("/home-link/iivow01/git/modules/work/26/4f24d38fa731c809caf0fcd2f85808/SRP254919.salmon.merged.gene_counts.top1000cov.tsv")

    GPROFILER2_GOST (
        input,
        background
    )
}


workflow test_gprofiler2_gost_backgroundlist {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true)
    ]

    GPROFILER2_GOST (
        input,
        []
    )
}
