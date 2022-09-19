#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRAFASTQ } from '../../../../subworkflows/nf-core/srafastq/main.nf'

workflow test_srafastq_single_end {
    input = Channel.of(
        [
            [ id:'test_single_end1', single_end:true ], // meta map
            'DRR000774'
        ],
        [
            [ id:'test_single_end2', single_end:true ], // meta map
            'DRR000775'
        ]
    )

    SRAFASTQ ( input )
}

workflow test_srafastq_paired_end {
    input = [
        [ id:'test_paired_end', single_end:false ], // meta map
        'SRR11140744'
    ]

    SRAFASTQ ( input )
}
