#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRODIGAL } from '../../../software/prodigal/main.nf' addParams( options: [:] )

workflow test_prodigal {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true) 
            ]

    PRODIGAL ( input , "gff")
}
