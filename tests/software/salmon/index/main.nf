#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../../software/salmon/index/main.nf' addParams( options: [:] )

workflow test_salmon_index {
    def genome_fasta     = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def transcript_fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_transcriptome.fasta", checkIfExists: true)
    SALMON_INDEX ( genome_fasta, transcript_fasta )
}
