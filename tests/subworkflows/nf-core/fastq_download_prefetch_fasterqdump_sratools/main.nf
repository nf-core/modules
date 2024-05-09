#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS } from '../../../../subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools/main.nf'

workflow test_fastq_download_prefetch_fasterqdump_sratools_single_end {
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

    FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS ( input, [] )
}

workflow test_fastq_download_prefetch_fasterqdump_sratools_paired_end {
    input = [
        [ id:'test_paired_end', single_end:false ], // meta map
        'SRR11140744'
    ]

    FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS ( input, [] )
}
