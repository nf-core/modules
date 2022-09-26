#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SELECTVARIANTS } from '../../../../modules/gatk4/selectvariants/main.nf'

// Basic parameters with uncompressed VCF input
workflow test_gatk4_selectvariants_vcf_input {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_idx'], checkIfExists: true)
    ]

    GATK4_SELECTVARIANTS ( input)
}

// Basic parameters with compressed VCF input
workflow test_gatk4_selectvariants_gz_input {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true)
    ]

    GATK4_SELECTVARIANTS ( input )
}
