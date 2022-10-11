#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_MARKDUPLICATES_PICARD } from '../../../../subworkflows/nf-core/bam_markduplicates_picard/main.nf'

workflow test_bam_markduplicates_picard {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BAM_MARKDUPLICATES_PICARD ( input, fasta )
}
