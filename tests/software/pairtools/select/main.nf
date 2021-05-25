#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process NO_PG {
    publishDir "${params.outdir}/nopg",
        mode: params.publish_dir_mode

    input:
      path sam
    output:
      path "${sam}.filtered", emit: fil
    script:
      """
      zcat < ${sam} | grep -v "@PG" > ${sam}.filtered
      """
}

include { PAIRTOOLS_SELECT } from '../../../../software/pairtools/select/main.nf' addParams( options: [args:"(pair_type == 'RU') or (pair_type == 'UR') or (pair_type == 'UU')", args2:'selected', args3:'rest'] )

workflow test_pairtools_select {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/open2c/pairtools/master/tests/data/mock.pairsam", checkIfExists: true) ]

    PAIRTOOLS_SELECT ( input )
    NO_PG(PAIRTOOLS_SELECT.out.sel.map{it[1]}.mix(PAIRTOOLS_SELECT.out.rest.map{it[1]}))
}
