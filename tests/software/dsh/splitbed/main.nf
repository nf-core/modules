#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSH_SPLITBED } from '../../../software/dsh/splitbed/main.nf' addParams( options: [suffix: '.', args: '--records 2'] )

workflow test_dsh_splitbed {

    def input = []
    input = [ [ id:'A' ], // meta map
              [ file("${launchDir}/tests/data/bed/A.bed.gz", checkIfExists: true) ] ]

    DSH_SPLITBED ( input )
}
