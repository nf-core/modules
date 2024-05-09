#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ_SOMATIC } from '../../../../../modules/nf-core/lofreq/somatic/main.nf'

workflow test_lofreq_somatic {

    input = [
        [ id:'test_', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    ]

    fasta = [
        [id:'test_genome'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:'test_genome_index'],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    LOFREQ_SOMATIC ( input, fasta, fai )
}
