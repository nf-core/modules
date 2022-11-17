#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONCOCNV } from '../../../../modules/nf-core/oncocnv/main.nf'

workflow test_oncocnv {
    
    normal = [
        "normal",
        [
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)
        ],
        [
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)
        ],
    ]
    tumor = [
        "tumor",
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_markduplicates_sorted_bam'], checkIfExists: true)
        ],
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true)
        ],
    ]

    ONCOCNV ( normal, tumor, file(params.test_data['homo_sapiens']['genome']['genome_21_annotated_bed'], checkIfExists: true), file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true) )
}
