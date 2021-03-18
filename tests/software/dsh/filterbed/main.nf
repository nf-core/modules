#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSH_FILTERBED } from '../../../../software/dsh/filterbed/main.nf' addParams( options: [suffix: '.filtered', args: '--range chr1:0-1000'] )

workflow test_dsh_filterbed {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true) ] ]

    DSH_FILTERBED ( input )
}
