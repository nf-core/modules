#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../../modules/salmon/index/main.nf'

workflow test_salmon_index {
    genome_fasta     = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    transcript_fasta = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)

    SALMON_INDEX ( genome_fasta, transcript_fasta )
}
