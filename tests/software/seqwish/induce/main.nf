#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQWISH_INDUCE } from '../../../../software/seqwish/induce/main.nf' addParams( options: [:] )

workflow test_seqwish_induce {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/genome/transcriptome.paf", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/genomics/sarscov2/genome/transcriptome.fasta", checkIfExists: true) ] 
            ]

    SEQWISH_INDUCE ( input )
}
