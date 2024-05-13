#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_DEEPVARIANT } from '../../../../../modules/nf-core/parabricks/deepvariant/main.nf'
include { PARABRICKS_DEEPVARIANT as PARABRICKS_DEEPVARIANT_GVCF } from '../../../../../modules/nf-core/parabricks/deepvariant/main.nf'

workflow test_parabricks_deepvariant {

    input = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        [], // bam not needed unless using intervals
        []
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    PARABRICKS_DEEPVARIANT ( input, fasta )
}

workflow test_parabricks_deepvariant_intervals {

    input = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]


    PARABRICKS_DEEPVARIANT ( input, fasta )
}

workflow test_parabricks_deepvariant_gvcf {

    input = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    PARABRICKS_DEEPVARIANT_GVCF ( input, fasta )
}