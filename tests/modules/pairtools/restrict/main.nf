#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_RESTRICT } from '../../../../modules/pairtools/restrict/main.nf'

workflow test_pairtools_restrict {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.4flip.pairs", checkIfExists: true) ]
    File dig  = new File("${workflow.workDir}/frag.bed")
    dig.write("chr1\t0\t50\r\nchr1\t50\t100\r\nchr2\t0\t50\r\nchr2\t50\t100\r\n!\t0\t1\r\n")
    frag  = file(dig)

    PAIRTOOLS_RESTRICT ( input, frag )
}
