#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADAM_MARKDUPLICATES } from '../../../../software/adam/markduplicates/main.nf' addParams( options: [:] )

workflow test_adam_markduplicates {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ] ]

    ADAM_MARKDUPLICATES ( input )
}
