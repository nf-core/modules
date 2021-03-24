#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../software/bwa/index/main.nf' addParams( options: [:] )

workflow test_bwa_index {
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BWA_INDEX ( fasta )
}
