#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_VARIANTS } from '../../../../modules/ivar/variants/main.nf'

workflow test_ivar_variants_no_gff_no_mpileup {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    gff   = []
    save_mpileup = false

    IVAR_VARIANTS ( input, fasta, fai, gff, save_mpileup )
}

workflow test_ivar_variants_no_gff_with_mpileup {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    gff   = []
    save_mpileup = true

    IVAR_VARIANTS ( input, fasta, fai, gff, save_mpileup )
}

workflow test_ivar_variants_with_gff_with_mpileup {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    gff   = file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true)
    save_mpileup = true

    IVAR_VARIANTS ( input, fasta, fai, gff, save_mpileup )
}
