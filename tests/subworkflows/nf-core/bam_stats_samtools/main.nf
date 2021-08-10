#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_stats_samtools/main' addParams( options: [:] )

workflow test_bam_stats_samtools_single_end {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
            ]

     BAM_STATS_SAMTOOLS ( input.join(SAMTOOLS_INDEX.out.bai, by: [0]) )
}

workflow test_bam_stats_samtools_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    BAM_STATS_SAMTOOLS ( input.join(SAMTOOLS_INDEX.out.bai, by: [0]) )
}
