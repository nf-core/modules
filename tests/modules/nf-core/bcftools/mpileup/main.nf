#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP } from '../../../../../modules/nf-core/bcftools/mpileup/main.nf'

workflow test_bcftools_mpileup {

    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        []
    ]

    fasta = [
        [ id:'sarscov2' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    ]

    save_mpileup = false

    BCFTOOLS_MPILEUP ( input, fasta, save_mpileup )
}

workflow test_bcftools_save_mpileup {

    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        []
    ]

    fasta = [
        [ id:'sarscov2' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    ]

    save_mpileup = true

    BCFTOOLS_MPILEUP ( input, fasta, save_mpileup )
}

workflow test_bcftools_mpileup_intervals {

    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]
    ]

    fasta = [
        [ id:'sarscov2' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    ]

    save_mpileup = false

    BCFTOOLS_MPILEUP ( input, fasta, save_mpileup )
}
