#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_RESTRICT } from '../../../../software/pairtools/restrict/main.nf' addParams( options: [:] )

workflow test_pairtools_restrict {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.4flip.pairs", checkIfExists: true) ]
    File dig  = new File("frag.bed")
    dig.write("chr1\t0\t50\nchr1\t50\t100\nchr2\t0\t50\nchr2\t50\t100\n!\t0\t1")
    frag  = file(dig)

    PAIRTOOLS_RESTRICT ( input, frag )
}
