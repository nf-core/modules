#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECLIP } from '../../../../modules/nf-core/pureclip/main.nf'

workflow test_pureclip_no_ctrl {
    input_bam = [ [ id:'test', single_end:false ], // meta map
              [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
              [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)  ]] //dummy file

    input_bai = [ [ id:'test', single_end:false ], // meta map
              [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ],
              [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true) ]] //dummy file

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    PURECLIP (
        input_bam,
        input_bai,
        fasta,
        false
        )
}

workflow test_pureclip_input_ctrl {
    input_bam = [ [ id:'test', single_end:false ], // meta map
              [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
              [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)  ]]

    input_bai = [ [ id:'test', single_end:false ], // meta map
              [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ],
              [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true) ]]

    genome = [file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]

    PURECLIP (
        input_bam,
        input_bai,
        genome,
        true
        )
}