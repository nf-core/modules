#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COLLECTSVEVIDENCE } from '../../../../../modules/nf-core/gatk4/collectsvevidence/main.nf'

workflow test_gatk4_collectsvevidence_cram {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_COLLECTSVEVIDENCE ( input, fasta, fasta_fai, dict )
}

workflow test_gatk4_collectsvevidence_bam {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        [],
        []
    ]

    fasta = []
    fasta_fai = []
    dict = []

    GATK4_COLLECTSVEVIDENCE ( input, fasta, fasta_fai, dict )
}

workflow test_gatk4_collectsvevidence_allele_count {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true)
    ]

    fasta = []
    fasta_fai = []
    dict = []

    GATK4_COLLECTSVEVIDENCE ( input, fasta, fasta_fai, dict )
}
