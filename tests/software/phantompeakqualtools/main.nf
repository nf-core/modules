#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHANTOMPEAKQUALTOOLS } from '../../../software/phantompeakqualtools/main.nf' addParams( options: [:] )

workflow test_phantompeakqualtools {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/bam/test.single_end.sorted.bam", checkIfExists: true) ] ]
    PHANTOMPEAKQUALTOOLS ( input )
}
