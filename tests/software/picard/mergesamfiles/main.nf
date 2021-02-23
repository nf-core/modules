#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MERGESAMFILES } from '../../../../software/picard/mergesamfiles/main.nf' addParams( options: [:] )

workflow test_picard_mergesamfiles {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3.bam", checkIfExists: true), ] ]

    PICARD_MERGESAMFILES ( input )
}
