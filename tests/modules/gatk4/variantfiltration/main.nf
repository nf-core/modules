#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

test_options = ['args': '--filter-name "test_filter" --filter-expression "MQ0 > 0"', 'suffix': '.filtered']
include { GATK4_VARIANTFILTRATION } from '../../../../modules/gatk4/variantfiltration/main.nf' addParams( options: test_options )

// Basic parameters with uncompressed VCF input
workflow test_gatk4_variantfiltration_vcf_input {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_idx'], checkIfExists: true) ]

    fasta        = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fastaIndex   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    fastaDict    = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_VARIANTFILTRATION ( input, fasta, fastaIndex, fastaDict )
}

// Basic parameters with compressed VCF input
workflow test_gatk4_variantfiltration_gz_input {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true) ]

    fasta        = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fastaIndex   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    fastaDict    = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_VARIANTFILTRATION ( input, fasta, fastaIndex, fastaDict )
}


