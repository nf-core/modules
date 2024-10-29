#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_MARKDUPLICATES_SAMTOOLS } from '../../../../subworkflows/nf-core/bam_markduplicates_samtools/main.nf'

workflow test_bam_markduplicates_samtools {

    input = [ [ id:'test' ], file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]

    fasta = [ [ id:'fasta' ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    BAM_MARKDUPLICATES_SAMTOOLS ( input, fasta )
}
