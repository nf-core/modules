#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GETPILEUPSUMMARIES } from '../../../../modules/gatk4/getpileupsummaries/main.nf' addParams( options: [:] )

workflow test_gatk4_getpileupsummaries_just_variants {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true) ,
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true) ]

    variants = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
    variants_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true)
    sites = []

    GATK4_GETPILEUPSUMMARIES ( input , variants , variants_tbi , sites )
}

workflow test_gatk4_getpileupsummaries_separate_sites {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true) ,
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true) ]

    variants = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
    variants_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true)
    sites = file( "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.interval_list" , checkIfExists: true)

    GATK4_GETPILEUPSUMMARIES ( input , variants , variants_tbi , sites )
}
