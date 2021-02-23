#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MARKDUPLICATES } from '../../../../software/picard/markduplicates/main.nf' addParams( options: [:] )

workflow test_picard_markduplicates_sorted_bam  {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]

    PICARD_MARKDUPLICATES ( input )
}

workflow test_picard_markduplicates_unsorted_bam  {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/sarscov2_paired_aln.bam", checkIfExists: true) ]

    PICARD_MARKDUPLICATES ( input )
}
