#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_STATS_SAMTOOLS } from '../../../../subworkflows/nf-core/bam_stats_samtools/main' addParams( options: [:] )

workflow test_bam_stats_samtools_single_end {
    input = [ [ id:'test', single_end:false ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]

    BAM_STATS_SAMTOOLS ( input )
}

workflow test_bam_stats_samtools_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    BAM_STATS_SAMTOOLS ( input )
}
