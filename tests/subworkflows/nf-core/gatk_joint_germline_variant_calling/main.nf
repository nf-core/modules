#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_JOINT_GERMLINE_VARIANT_CALLING } from '../../../../subworkflows/nf-core/gatk_joint_germline_variant_calling/main'

workflow test_gatk_joint_germline_variant_calling {
    input = [
             [ [ id:'test' ], // meta map
               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
             ],
             [ [ id:'test2' ], // meta map
               file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
               file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
             ]
            ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    sites = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    sites_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)
    intervals = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    joint_id = "test_joint"
    GATK_JOINT_GERMLINE_VARIANT_CALLING ( input, fasta, fai, dict, sites, sites_tbi, intervals, joint_id )
}
