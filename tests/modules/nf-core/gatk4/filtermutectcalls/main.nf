#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FILTERMUTECTCALLS } from '../../../../modules/gatk4/filtermutectcalls/main.nf'

workflow test_gatk4_filtermutectcalls_base {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_stats'], checkIfExists: true),
        [],
        [],
        [],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    GATK4_FILTERMUTECTCALLS ( input, fasta, fai, dict )
}

workflow test_gatk4_filtermutectcalls_with_files {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_stats'], checkIfExists: true),
        [ file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_artifact_prior_tar_gz'], checkIfExists: true) ],
        [ file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_segmentation_table'], checkIfExists: true) ],
        [ file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_contamination_table'], checkIfExists: true) ],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    GATK4_FILTERMUTECTCALLS ( input, fasta, fai, dict )
}

workflow test_gatk4_filtermutectcalls_use_val {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_stats'], checkIfExists: true),
        [ file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_artifact_prior_tar_gz'], checkIfExists: true) ],
        [ file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_segmentation_table'], checkIfExists: true) ],
        [],
        '20.0'
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    GATK4_FILTERMUTECTCALLS ( input, fasta, fai, dict )
}
