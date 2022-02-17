#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRELKA_GERMLINE } from '../../../../modules/strelka/germline/main.nf'

workflow test_strelka_germline {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    STRELKA_GERMLINE ( input, fasta, fai)
}

workflow test_strelka_germline_target_bed {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed_gz_tbi'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    STRELKA_GERMLINE ( input, fasta, fai)
}

