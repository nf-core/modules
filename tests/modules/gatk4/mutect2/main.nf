#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MUTECT2 } from '../../../../modules/gatk4/mutect2/main.nf' addParams( options: [:] )

workflow test_gatk4_mutect2_tumor_normal_pair {

input = [ [ id:'test' , run_pon:false, run_single:false ], // meta map
           [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true) , file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
           [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true) , file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
           ["testN"]
           ]


fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
germline_resource_idx = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true)
panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true)
panel_of_normals_idx = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'], checkIfExists: true)


    GATK4_MUTECT2 ( input , fasta , fastaidx , dict , germline_resource, germline_resource_idx , panel_of_normals , panel_of_normals_idx )
}

workflow test_gatk4_mutect2_tumor_single {

input = [ [ id:'test' , run_pon:false, run_single:true ], // meta map
           [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
           [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
           []
           ]


fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
germline_resource = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
germline_resource_idx = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true)
panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true)
panel_of_normals_idx = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'], checkIfExists: true)


    GATK4_MUTECT2 ( input , fasta , fastaidx , dict , germline_resource, germline_resource_idx , panel_of_normals , panel_of_normals_idx )
}

workflow test_gatk4_mutect2_generate_pon {

input = [ [ id:'test' , run_pon:true, run_single:false ], // meta map
           [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
           [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
           []
           ]


fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
germline_resource = []
germline_resource_idx = []
panel_of_normals = []
panel_of_normals_idx = []


    GATK4_MUTECT2 ( input , fasta , fastaidx , dict , germline_resource, germline_resource_idx , panel_of_normals , panel_of_normals_idx )
}
