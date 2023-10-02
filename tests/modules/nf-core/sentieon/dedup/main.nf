#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_DEDUP as SENTIEON_DEDUP_MARK   } from '../../../../../modules/nf-core/sentieon/dedup/main.nf'
include { SENTIEON_DEDUP as SENTIEON_DEDUP_REMOVE } from '../../../../../modules/nf-core/sentieon/dedup/main.nf'

workflow test_dedup_mark_duplicate_reads {

    fasta_file = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai_file = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    bam_ch = [
        [ id: 'test' ],
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]

    SENTIEON_DEDUP_MARK ( bam_ch, [[:], fasta_file], [[:], fasta_fai_file] )
}

workflow test_dedup_remove_duplicate_reads {

    fasta_file = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai_file = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    bam_ch = [
        [ id: 'test' ],
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]

    SENTIEON_DEDUP_REMOVE ( bam_ch, [[:], fasta_file], [[:], fasta_fai_file] )
}
