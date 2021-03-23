#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENZA_WIGGLE } from '/home/AD/rbhuller/modules/software/sequenza/wiggle/main.nf' addParams( options: [ 'args': '-w 50' ] )

workflow test_sequenza_wiggle {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]

    SEQUENZA_WIGGLE ( input )
}
