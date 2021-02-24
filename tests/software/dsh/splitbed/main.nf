#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSH_SPLITBED } from '../../../../software/dsh/splitbed/main.nf' addParams( options: [suffix: '.', args: '--records 2'] )

workflow test_dsh_splitbed {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true) ] ]

    DSH_SPLITBED ( input )
}
