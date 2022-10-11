#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_MARKDUPLICATES_PICARD; BAM_MARKDUPLICATES_PICARD as BAM_MARKDUPLICATES_PICARD_CRAM } from '../../../../subworkflows/nf-core/bam_markduplicates_picard/main.nf'

workflow test_bam_markduplicates_picard {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BAM_MARKDUPLICATES_PICARD ( input, fasta )
}

workflow test_bam_markduplicates_picard_cram {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BAM_MARKDUPLICATES_PICARD_CRAM ( input, fasta )
}
