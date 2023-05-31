#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../../../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'

workflow test_bam_dedup_stats_samtools_umitools {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    get_output_stats = false

    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS ( input, get_output_stats )
}
