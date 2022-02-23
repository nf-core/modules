#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONTROLFREEC_SOMATIC } from '../../../modules/controlfreec/main.nf'
include { SAMTOOLS_MPILEUP } from '../../../modules/samtools/mpileup/main.nf'

workflow test_controlfreec {

    input = [
        [ id:'test', single_end:false, sex:'XX' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        [], [], [], []
    ]




    CONTROLFREEC_SOMATIC (input, [], [], [], [], [], [], [], [], [])
}
