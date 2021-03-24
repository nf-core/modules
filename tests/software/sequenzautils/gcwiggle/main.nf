#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENZAUTILS_GCWIGGLE } from '../../../../software/sequenzautils/gcwiggle/main.nf' addParams( options: [ 'args': '-w 50' ] )

workflow test_sequenzautils_gcwiggle {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]

    SEQUENZAUTILS_GCWIGGLE ( input )
}
