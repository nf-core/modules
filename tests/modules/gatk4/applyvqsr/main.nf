#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_APPLYVQSR } from '../../../../modules/gatk4/applyvqsr/main.nf'

workflow test_gatk4_applyvqsr {
    input = [ [ id:'test'], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_recal'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_recal_idx'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_tranches'], checkIfExists: true)
            ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    allelespecific = false
    truthsensitivity = '99.0'
    mode = 'SNP'

    GATK4_APPLYVQSR ( input, fasta, fai, dict, allelespecific, truthsensitivity, mode )
}

workflow test_gatk4_applyvqsr_allele_specific {
    input = [ [ id:'test'], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_allele_specific_recal'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_allele_specific_recal_idx'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_allele_specific_tranches'], checkIfExists: true)
            ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    allelespecific = true
    truthsensitivity = '99.0'
    mode = 'SNP'

    GATK4_APPLYVQSR ( input, fasta, fai, dict, allelespecific, truthsensitivity, mode )
}
