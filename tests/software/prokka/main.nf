#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROKKA } from '../../../software/prokka/main.nf' addParams( options: [:] )

workflow test_prokka {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
            ]
    
    PROKKA ( input, [], [] )
}
