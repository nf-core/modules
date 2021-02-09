#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TOOL_SUBTOOL } from '../../../../software/TOOL/SUBTOOL/main.nf' addParams( options: [:] )

workflow test_tool_subtool {
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]

    TOOL_SUBTOOL ( input )
}