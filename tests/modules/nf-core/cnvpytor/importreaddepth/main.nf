#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVPYTOR_IMPORTREADDEPTH } from '../../../../modules/cnvpytor/importreaddepth/main.nf'


workflow test_cnvpytor_importreaddepth {

    input = [
        [ id: 'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    CNVPYTOR_IMPORTREADDEPTH (input, [], [])
}

workflow test_cnvpytor_importreaddepth_cram {

    input = [
        [ id: 'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    CNVPYTOR_IMPORTREADDEPTH (input, fasta, fai)
}
