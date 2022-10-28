#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTA_BINNING_CONCOCT } from '../../../../subworkflows/nf-core/fasta_binning_concoct/main.nf'

workflow test_fasta_binning_concoct {

    ch_fasta = Channel.of(
            [[ id: 'test' ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)],
            [[ id: 'test2'], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
        ).dump(tag: "fasta")
    ch_bam   = Channel.of(
            [[ id: 'test' ], file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)],
            [[ id: 'test2'], file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]
        ).dump(tag: "bam")

    FASTA_BINNING_CONCOCT ( ch_fasta, ch_bam )
}
