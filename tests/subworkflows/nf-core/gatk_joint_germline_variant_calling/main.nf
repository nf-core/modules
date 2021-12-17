#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_JOINT_GERMLINE_VARIANT_CALLING } from '../../../../subworkflows/nf-core/gatk_joint_germline_variant_calling/main'

workflow test_gatk_joint_germline_variant_calling_skip_haplotc {
    input            = [
                        [[ id:'test' ], // meta map
                         file(params.test_data['homo_sapiens']['illumina']['test_g_vcf_gz']                         ,  checkIfExists: true),
                         file(params.test_data['homo_sapiens']['illumina']['test_g_vcf_gz_tbi']                     ,  checkIfExists: true)
                        ],
                        [[ id:'test2' ], // meta map
                         file(params.test_data['homo_sapiens']['illumina']['test2_g_vcf_gz']                        , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['illumina']['test2_g_vcf_gz_tbi']                    , checkIfExists: true)
                        ]
                       ]
    run_haplotc      = false
    run_vqsr         = true
    fasta            = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta']                           , checkIfExists: true)
    fai              = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai']                       , checkIfExists: true)
    dict             = file(params.test_data['homo_sapiens']['genome']['genome_21_dict']                            , checkIfExists: true)
    sites            = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz']                  , checkIfExists: true)
    sites_tbi        = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi']              , checkIfExists: true)
    intervals        = file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list']                   , checkIfExists: true)
    joint_id         = "test_joint"
    allelespecific   = false
    resources        = [
                        [
                         file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz']               , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz']        , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz']    , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz']                , checkIfExists: true)
                        ],
                        [
                         file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz_tbi']           , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz_tbi']    , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz_tbi'], checkIfExists: true),
                         file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi']            , checkIfExists: true)
                        ],
                        [
                         'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
                         'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
                         '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.hg38.vcf.gz',
                         'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf.gz'
                        ]
                       ]
    annotation       = ['QD', 'FS', 'SOR']
    mode             = 'SNP'
    create_rscript   = false
    truthsensitivity = '99.0'
    GATK_JOINT_GERMLINE_VARIANT_CALLING ( input, run_haplotc, run_vqsr, fasta, fai, dict, sites, sites_tbi, intervals, joint_id, allelespecific, resources, annotation, mode, create_rscript, truthsensitivity )
}

workflow test_gatk_joint_germline_variant_calling_skip_vqsr {
    input            = [
                        [[ id:'test' ], // meta map
                         file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam']     , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'] , checkIfExists: true)
                        ],
                        [[ id:'test2' ], // meta map
                         file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam']    , checkIfExists: true),
                         file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
                        ]
                       ]
    run_haplotc      = true
    run_vqsr         = false
    fasta            = file(params.test_data['homo_sapiens']['genome']['genome_fasta']                       , checkIfExists: true)
    fai              = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai']                   , checkIfExists: true)
    dict             = file(params.test_data['homo_sapiens']['genome']['genome_dict']                        , checkIfExists: true)
    sites            = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz']              , checkIfExists: true)
    sites_tbi        = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi']          , checkIfExists: true)
    intervals        = file(params.test_data['homo_sapiens']['genome']['genome_bed']                         , checkIfExists: true)
    joint_id         = "test_joint"
    allelespecific   = []
    resources        = []
    annotation       = []
    mode             = []
    create_rscript   = []
    truthsensitivity = []
    GATK_JOINT_GERMLINE_VARIANT_CALLING ( input, run_haplotc, run_vqsr, fasta, fai, dict, sites, sites_tbi, intervals, joint_id, allelespecific, resources, annotation, mode, create_rscript, truthsensitivity )
}
