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

include { PAIRTOOLS_SORT } from '../../../../software/pairtools/sort/main.nf' addParams( options: [:] )

workflow test_pairtools_sort {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.pairsam", checkIfExists: true) ]

    PAIRTOOLS_SORT ( input ).sorted.map{it[1]} | NO_PG
}
