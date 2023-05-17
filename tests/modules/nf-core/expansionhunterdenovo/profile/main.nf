#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXPANSIONHUNTERDENOVO_PROFILE } from '../../../../../modules/nf-core/expansionhunterdenovo/profile/main.nf'

workflow test_expansionhunterdenovo_profile {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)


    EXPANSIONHUNTERDENOVO_PROFILE (
        input,
        fasta,
        fasta_fai
    )
}

workflow test_expansionhunterdenovo_profile_cram {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)


    EXPANSIONHUNTERDENOVO_PROFILE (
        input,
        fasta,
        fasta_fai
    )
}

