#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process NO_PG {
    publishDir "${params.outdir}/nopg",
        mode: params.publish_dir_mode

    input:
      path sam
    output:
      path "filtered.file", emit: fil
    script:
      """
      zcat < ${sam} | grep -v "@PG" > filtered.file
      """
}

include { PAIRTOOLS_RESTRICT } from '../../../../software/pairtools/restrict/main.nf' addParams( options: [:] )

workflow test_pairtools_restrict {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.4flip.pairs", checkIfExists: true) ]
    File dig  = new File("frag.bed")
    dig.write("chr1\t0\t50\nchr1\t50\t100\nchr2\t0\t50\nchr2\t50\t100\n!\t0\t1")
    frag  = file(dig)

    PAIRTOOLS_RESTRICT ( input, frag ).restrict.map{it[1]} | NO_PG
}
