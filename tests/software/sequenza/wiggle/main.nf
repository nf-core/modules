#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENZA_WIGGLE } from '../../../software/sequenza/wiggle/main.nf' addParams( options: [ 'args': '-w 50' ] )

workflow test_sequenza {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]

    SEQUENZA_WIGGLE ( input )
}
