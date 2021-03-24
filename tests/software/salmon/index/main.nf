#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../../software/salmon/index/main.nf' addParams( options: [:] )

workflow test_salmon_index {
    genome        = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    transcriptome = file("${launchDir}/tests/data/genomics/sarscov2/genome/transcriptome.fasta", checkIfExists: true)

    SALMON_INDEX ( genome, transcriptome )
}
