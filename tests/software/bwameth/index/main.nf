#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMETH_INDEX } from '../../../../software/bwameth/index/main.nf' addParams( options: [:] )

workflow test_bwameth_index {

    BWAMETH_INDEX ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
}
