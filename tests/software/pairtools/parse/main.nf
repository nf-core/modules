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

include { PAIRTOOLS_PARSE } from '../../../../software/pairtools/parse/main.nf' addParams( options: [:] )

workflow test_pairtools_parse {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.sam", checkIfExists: true) ]
    sizes = file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.chrom.sizes", checkIfExists:true)

    PAIRTOOLS_PARSE ( input, sizes ).pairsam.map{it[1]} | NO_PG
}
