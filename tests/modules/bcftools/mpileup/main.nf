#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP } from '../../../../modules/bcftools/mpileup/main.nf'

workflow test_bcftools_mpileup {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]]
    fasta = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    BCFTOOLS_MPILEUP ( input, fasta )
}
