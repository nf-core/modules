#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTO_INDEX } from '../../../../software/kallisto/index/main.nf' addParams( options: [:] )

workflow test_kallisto_index {

    def input = []
    input = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    KALLISTO_INDEX ( input )
}
