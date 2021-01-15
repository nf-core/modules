#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRINGTIE as STRINGTIE_FORWARD } from '../../../software/stringtie/main.nf'  addParams( options: [ publish_dir:'test_stringtie_forward' ] )
include { STRINGTIE as STRINGTIE_REVERSE } from '../../../software/stringtie/main.nf'  addParams( options: [ publish_dir:'test_stringtie_reverse' ] )

/*
 * Test with forward strandedness
 */
workflow test_stringtie_forward {
    def input = []
    input = [ [ id:'test', strandedness:'forward' ], // meta map
              [ file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ] ]

    STRINGTIE_FORWARD (
        input,
        file("${launchDir}/tests/data/gff/a.gtf", checkIfExists: true)
    )
}

/*
 * Test with reverse strandedness
 */
workflow test_stringtie_reverse {
    def input = []
    input = [ [ id:'test', strandedness:'reverse' ], // meta map
              [ file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ] ]

    STRINGTIE_REVERSE (
        input,
        file("${launchDir}/tests/data/gff/a.gtf", checkIfExists: true)
    )
}
