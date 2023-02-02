#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLTOOLS_PILEUP } from '../../../../../modules/nf-core/cooltools/pileup/main.nf'

workflow test_cooltools_pileup {

    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/open2c/cooltools/master/tests/data/CN.mm9.1000kb.cool', checkIfExists: true),
        [:] // resolution if any
    ]
    File bed  = new File("${workflow.workDir}/test.bed")
    bed.write("chr1\t102000000\t107000000\r\nchr1\t108000000\t113000000\r\n")
    frag  = file(bed)

    COOLTOOLS_PILEUP ( input, frag )
}
