#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MUTECT2 } from '../../../../modules/gatk4/mutect2/main.nf' addParams( options: [:] )
// used to run with the mitochondria mode setting as this increases sensitivity, allowing for some tumor_normal variants to be detected while the old test data is still in use, will be removed when new test data for sarek is available.
include { GATK4_MUTECT2 as GATK4_TEMPFIX_MUTECT2 } from '../../../../modules/gatk4/mutect2/main.nf' addParams( options: [args: '--mitochondria-mode'] )

workflow test_gatk4_mutect2_tumor_normal_pair {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
              ["testN"]
            ]
    run_single = false
    run_pon = false
    run_mito = false
    interval_label = []
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'], checkIfExists: true)

    GATK4_TEMPFIX_MUTECT2 ( input, run_single, run_pon, run_mito, interval_label, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

workflow test_gatk4_mutect2_tumor_single {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
              []
            ]
    run_single = true
    run_pon = false
    run_mito = false
    interval_label = []
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'], checkIfExists: true)

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, interval_label, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

workflow test_gatk4_mutect2_cram_input {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true)],
              []
            ]
    run_single = true
    run_pon = false
    run_mito = false
    interval_label = []
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'], checkIfExists: true)

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, interval_label, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

workflow test_gatk4_mutect2_generate_pon {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
              []
            ]
    run_single = false
    run_pon = true
    run_mito = false
    interval_label = []
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    germline_resource = []
    germline_resource_tbi = []
    panel_of_normals = []
    panel_of_normals_tbi = []

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, interval_label, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

// mitochondria mode would ideally have some mitochondria test data, but since the mitochondria settings only increase detection sensitivity, we can use the chr22 data as a stand in as it is already a small dataset, the extra variants detected compared to generate_pon shows the mode is working.
workflow test_gatk4_mutect2_mitochondria {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
              []
            ]
    run_single = false
    run_pon = false
    run_mito = true
    interval_label = 'chr22'
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    germline_resource = []
    germline_resource_tbi = []
    panel_of_normals = []
    panel_of_normals_tbi = []

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, interval_label, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}
