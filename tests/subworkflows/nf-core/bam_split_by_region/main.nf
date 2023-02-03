#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../../modules/nf-core/samtools/index/main.nf'
include { BAM_SPLIT_BY_REGION } from '../../../../subworkflows/nf-core/bam_split_by_region/main.nf'

workflow test_bam_split_by_region {

    input = [
        [ [ id: 'test', single_end: false ],
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'],     checkIfExists: true),
        ],
        [ [ id: 'test2', single_end: false ],
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'],     checkIfExists: true),
        ]
    ]

    regions = [
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'],     checkIfExists: true)
    ]

    SAMTOOLS_INDEX(Channel.fromList(input))

    ch_bam_split_input=Channel.fromList(input)
            .join(SAMTOOLS_INDEX.out.bai)
            .combine(Channel.of(regions))

    BAM_SPLIT_BY_REGION ( ch_bam_split_input )
}
