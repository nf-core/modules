#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GETPILEUPSUMMARIES } from '../../../../modules/gatk4/getpileupsummaries/main.nf'

workflow test_gatk4_getpileupsummaries_just_variants {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true) ,
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        []
        ]

    variants = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    variants_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    fasta = []
    fai = []
    dict = []
    GATK4_GETPILEUPSUMMARIES ( input , fasta, fai, dict, variants , variants_tbi )
}

workflow test_gatk4_getpileupsummaries_separate_sites {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true) ,
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true) ]

    variants = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    variants_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    fasta = []
    fai = []
    dict = []
    GATK4_GETPILEUPSUMMARIES ( input , fasta, fai, dict, variants , variants_tbi)
}

workflow test_gatk4_getpileupsummaries_separate_sites_cram {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true) ,
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true)
        ]

    variants = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    variants_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    GATK4_GETPILEUPSUMMARIES ( input , fasta, fai, dict, variants , variants_tbi)
}
