#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEACR_CALLPEAK } from '../../../../modules/seacr/callpeak/main.nf' addParams( options: [ args:'norm stringent' ] )

workflow test_seacr_callpeak {
    input = [
        [ id:'test_1'],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/bedgraph/K27me3_1_to_chr20.bedgraph", checkIfExists: true),
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/bedgraph/IgG_1_to_chr20.bedgraph", checkIfExists: true)
    ]

    SEACR_CALLPEAK ( input )
}
