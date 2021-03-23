#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROKKA } from '../../../software/prokka/main.nf' addParams( options: [:] )

workflow test_prokka {
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true), ]
    PROKKA ( input, [], [] )
}
