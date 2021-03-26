#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MERGE } from '../../../../software/samtools/merge/main.nf' addParams( options: [:] )

workflow test_samtools_merge {
    input = [ [ id: 'test' ], // meta map
               [ file(params.test_data['sarscov2']['illumina']['test_methylated_paired_end_sorted_bam'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)]
               ]

    SAMTOOLS_MERGE ( input )
}
