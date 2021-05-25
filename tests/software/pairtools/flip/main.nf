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

include { PAIRTOOLS_FLIP } from '../../../../software/pairtools/flip/main.nf' addParams( options: [:] )

workflow test_pairtools_flip {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.4flip.pairs", checkIfExists: true) ]
    sizes = file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.chrom.sizes", checkIfExists:true)

    PAIRTOOLS_FLIP ( input, sizes ).flip.map{it[1]} | NO_PG
}
