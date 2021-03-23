#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTO_INDEX } from '../../../../software/kallisto/index/main.nf' addParams( options: [:] )

workflow test_kallisto_index {
    
    def genome_fasta = file("${launchDir}/tests/data/genomics/sarscov2/bam/test_single_end.bam", checkIfExists: true)
    def transcript_fasta = file("${launchDir}/tests/data/genomics/sarscov2/bam/test_single_end.bam", checkIfExists: true)

    KALLISTO_INDEX ( genome_fasta, transcript_fasta )
}
