#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_GVCFTYPER } from '../../../../../modules/nf-core/sentieon/gvcftyper/main.nf'
include { UNTAR              } from '../../../../../modules/nf-core/untar/main.nf'

// Basic parameters with uncompressed VCF input
workflow test_sentieon_gvcftyper_vcf_input {

    input = [ [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_idx'], checkIfExists: true),
            []
        ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    SENTIEON_GVCFTYPER( input, fasta, fai, [], [])
}

// Basic parameters with compressed VCF input
workflow test_sentieon_gvcftyper_gz_input {

    input = [ [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
            []
        ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    SENTIEON_GVCFTYPER( input, fasta, fai, [], [])
}

// Basic parameters + optional dbSNP
workflow test_sentieon_gvcftyper_gz_input_dbsnp {

    input = [ [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
            []
        ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dbsnp     = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    dbsnp_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    SENTIEON_GVCFTYPER( input, fasta, fai, dbsnp, dbsnp_tbi)
}

// Basic parameters + optional intervals
workflow test_sentieon_gvcftyper_gz_input_intervals {

    input = [ [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)            ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    SENTIEON_GVCFTYPER( input, fasta, fai, [], [])
}

// Basic parameters + optional dbSNP + optional intervals
workflow test_sentieon_gvcftyper_gz_input_dbsnp_intervals {

    input = [ [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
        ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dbsnp     = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    dbsnp_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    SENTIEON_GVCFTYPER( input, fasta, fai, dbsnp, dbsnp_tbi )
}
