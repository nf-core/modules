#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_VARIANTFILTRATION } from '../../../../modules/gatk4/variantfiltration/main.nf'

// Basic parameters with uncompressed VCF input
workflow test_gatk4_variantfiltration_vcf_input {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_idx'], checkIfExists: true)
    ]

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_index = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    fasta_dict  = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_VARIANTFILTRATION ( input, fasta, fasta_index, fasta_dict )
}

// Basic parameters with compressed VCF input
workflow test_gatk4_variantfiltration_gz_input {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true)
    ]

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_index = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    fasta_dict  = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_VARIANTFILTRATION ( input, fasta, fasta_index, fasta_dict )
}


