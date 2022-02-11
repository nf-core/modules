#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSORPRO_MSI_SOMATIC } from '../../../../modules/msisensorpro/msi_somatic/main.nf'
include { MSISENSORPRO_SCAN } from '../../../../modules/msisensorpro/scan/main.nf'

workflow test_msisensorpro_msi {

    scan_in = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    println scan_in

    MSISENSORPRO_SCAN ( scan_in )

    input = [// meta map
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        []
    ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)

    MSISENSORPRO_SCAN.out.list.map{meta, list -> [list]}.set{list}
    MSISENSORPRO_MSI_SOMATIC(input, fasta, list)


}
