#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FILTERVARIANTTRANCHES } from '../../../../modules/gatk4/filtervarianttranches/main.nf'

workflow test_gatk4_filtervarianttranches {

    input = [ [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
                []
            ]

    resources = [
                    file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz'], checkIfExists: true),
                    //file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz'], checkIfExists: true),
                    //file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz'], checkIfExists: true),
                    //file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
                ]
    resources_index  =[              file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz_tbi'], checkIfExists: true),
                    //file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz_tbi'], checkIfExists: true),
                    //file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz_tbi'], checkIfExists: true),
                    //file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)
                ]
                // , [
                //     'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
                //     'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
                //     '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.hg38.vcf.gz',
                //     'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf.gz'
                // ]]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    GATK4_FILTERVARIANTTRANCHES ( input, resources, resources_index, fasta, fai, dict)
}
