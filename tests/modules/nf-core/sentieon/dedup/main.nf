#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_DEDUP } from '../../../../../modules/nf-core/sentieon/dedup/main.nf'

workflow test_sentieon_dedup {

    fasta_file = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai_file = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    bam_ch = [
        [ id: 'test' ],
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]

    SENTIEON_DEDUP ( bam_ch, fasta_file, fasta_fai_file )
}
