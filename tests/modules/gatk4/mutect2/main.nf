#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MUTECT2 } from '../../../../modules/gatk4/mutect2/main.nf'

workflow test_gatk4_mutect2_tumor_normal_pair {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)
                ],
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)
                ],
              [],
              ["normal"]
            ]
    run_single = false
    run_pon = false
    run_mito = false
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

workflow test_gatk4_mutect2_tumor_single {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
              [],
              []
            ]
    run_single = true
    run_pon = false
    run_mito = false
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

workflow test_gatk4_mutect2_cram_input {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true)],
              [],
              []
            ]
    run_single = true
    run_pon = false
    run_mito = false
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

workflow test_gatk4_mutect2_generate_pon {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
              [],
              []
            ]
    run_single = false
    run_pon = true
    run_mito = false
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    germline_resource = []
    germline_resource_tbi = []
    panel_of_normals = []
    panel_of_normals_tbi = []

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}

// mitochondria mode would ideally have some mitochondria test data, but since the mitochondria settings only increase detection sensitivity, we can use the chr22 data as a stand in as it is already a small dataset, the extra variants detected compared to generate_pon shows the mode is working.
workflow test_gatk4_mutect2_mitochondria {
    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['mitochon_standin_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['mitochon_standin_recalibrated_sorted_bam_bai'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)],
              []
            ]
    run_single = false
    run_pon = false
    run_mito = true
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    germline_resource = []
    germline_resource_tbi = []
    panel_of_normals = []
    panel_of_normals_tbi = []

    GATK4_MUTECT2 ( input, run_single, run_pon, run_mito, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
}
