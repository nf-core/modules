#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_CONSENSUS } from '../../../../modules/ivar/consensus/main.nf'

workflow test_ivar_consensus {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    save_mpileup = false

    IVAR_CONSENSUS ( input, fasta, save_mpileup)
}

workflow test_ivar_consensus_mpileup {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    save_mpileup = true

    IVAR_CONSENSUS ( input, fasta, save_mpileup)
}
