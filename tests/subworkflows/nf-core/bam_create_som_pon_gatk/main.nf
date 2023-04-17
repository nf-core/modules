#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_CREATE_SOM_PON_GATK } from '../../../../subworkflows/nf-core/bam_create_som_pon_gatk/main'

workflow bam_create_som_pon_gatk {
    ch_mutect2_in = Channel.of(
        [
            [ id:'test1' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
            []
        ],
        [
            [ id:'test2' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
            []
        ]
    )
    fasta           = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true))
    fai             = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true))
    dict            = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true))
    interval_file   = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true))
    pon_norm        = "test_panel"

    BAM_CREATE_SOM_PON_GATK ( ch_mutect2_in, fasta, fai, dict, pon_norm, interval_file )
}
