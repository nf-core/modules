#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../../software/salmon/index/main.nf' addParams( options: [:] )

workflow test_salmon_index {
    def genome_fasta     = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    def transcript_fasta = file(params.test_data['sarscov2']['genome']['transcriptome.fasta'], checkIfExists: true)
    SALMON_INDEX ( genome_fasta, transcript_fasta )
}
