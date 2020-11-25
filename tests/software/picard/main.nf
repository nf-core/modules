#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MERGESAMFILES } from '../../../software/picard/mergesamfiles/main.nf' addParams( options: [:] )

workflow test_picard_mergesamfiles {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true),
                file("${launchDir}/tests/data/bam/test.paired_end.name.sorted.bam", checkIfExists: true), ] ]

    PICARD_MERGESAMFILES ( input )
}
