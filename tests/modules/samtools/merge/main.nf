#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MERGE } from '../../../../modules/samtools/merge/main.nf'

workflow test_samtools_merge {
    input = [ [ id: 'test' ], // meta map
               [ file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)]
               ]

    SAMTOOLS_MERGE ( input, [] )
}

workflow test_samtools_merge_cram {
    input = [ [ id: 'test' ], // meta map
               [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
                 file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
                 ]
               ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    SAMTOOLS_MERGE ( input, fasta )
}
