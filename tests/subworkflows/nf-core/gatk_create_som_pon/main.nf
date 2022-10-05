#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_CREATE_SOM_PON } from '../../../../subworkflows/nf-core/gatk_create_som_pon/main' addParams( [:]  )

workflow test_gatk_create_som_pon {
    ch_mutect2_in = [
                    [[ id:'test1' ], // meta map
                    [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
                    [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
                    [] ],
                    [[ id:'test2' ], // meta map
                    [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
                    [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
                    [] ]
                    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    pon_name = "test_panel"
    interval_file = file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)

    GATK_CREATE_SOM_PON ( ch_mutect2_in, fasta, fai, dict, pon_name, interval_file )
}
