#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../../modules/nf-core/samtools/index/main.nf'
include { BAM_SPLIT_BY_CHROM } from '../../../../subworkflows/nf-core/bam_split_by_chrom/main.nf'

workflow test_bam_split_by_chrom {

    input = [ [ id: 'test', single_end: false ],
        file(params.test_data['homo_sapiens']['illumina']['test3_single_end_markduplicates_sorted_bam'],     checkIfExists: true)
    ]

    SAMTOOLS_INDEX(input)

    ch_bam_split_input=Channel.of(input)
            .join(SAMTOOLS_INDEX.out.bai)

    BAM_SPLIT_BY_CHROM ( ch_bam_split_input )
}
