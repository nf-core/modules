#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_READWRITER as SENTIEON_READWRITER_BAM  } from '../../../../../modules/nf-core/sentieon/readwriter/main.nf'
include { SENTIEON_READWRITER as SENTIEON_READWRITER_CRAM } from '../../../../../modules/nf-core/sentieon/readwriter/main.nf'

workflow test_readwriter_bam {

    bam = [
        [ id: 'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    SENTIEON_READWRITER_BAM ( bam, [[:],[]], [[:],[]] )
}

workflow test_readwriter_cram {

    cram = [
        [ id: 'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]

    fasta = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]).collect()
    fai   = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]).collect()

    SENTIEON_READWRITER_CRAM ( cram, fasta, fai )
}
