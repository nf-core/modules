#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENZAUTILS_GCWIGGLE } from '../../../../software/sequenzautils/gcwiggle/main.nf' addParams( options: [ 'args': '-w 50' ] )

workflow test_sequenzautils_gcwiggle {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true) 
            ]

    SEQUENZAUTILS_GCWIGGLE ( input )
}
