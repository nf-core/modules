#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARAGRAPH_IDXDEPTH } from '../../../../../modules/nf-core/paragraph/idxdepth/main.nf'

workflow test_paragraph_idxdepth_bam {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta = [
        [ id:'fasta' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'fasta_fai' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    PARAGRAPH_IDXDEPTH (
        input,
        fasta,
        fasta_fai
    )
}

workflow test_paragraph_idxdepth_cram {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]

    fasta = [
        [ id:'fasta' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'fasta_fai' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    PARAGRAPH_IDXDEPTH (
        input,
        fasta,
        fasta_fai
    )
}
