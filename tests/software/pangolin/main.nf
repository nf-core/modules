#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PANGOLIN } from '../../../software/pangolin/main.nf' addParams( options: [:] )

workflow test_pangolin {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ] ]

    PANGOLIN ( input )
}
