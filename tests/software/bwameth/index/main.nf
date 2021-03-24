#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMETH_INDEX } from '../../../../software/bwameth/index/main.nf' addParams( options: [:] )

workflow test_bwameth_index {
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BWAMETH_INDEX ( fasta )
}
