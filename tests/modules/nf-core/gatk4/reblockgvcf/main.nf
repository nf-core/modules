#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_REBLOCKGVCF } from '../../../../../modules/nf-core/gatk4/reblockgvcf/main.nf'

workflow test_gatk4_reblockgvcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        []
    ]

    fasta       = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_index = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict        = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_REBLOCKGVCF ( input, fasta, fasta_index, dict, [], [] )
}

workflow test_gatk4_reblockgvcf_intervals {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    fasta       = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_index = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict        = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_REBLOCKGVCF ( input, fasta, fasta_index, dict, [], [] )
}

workflow test_gatk4_reblockgvcf_dbsnp {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
        []
    ]

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_index = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict        = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    dbsnp       = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    dbsnp_tbi   = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    GATK4_REBLOCKGVCF ( input, fasta, fasta_index, dict, dbsnp, dbsnp_tbi )
}