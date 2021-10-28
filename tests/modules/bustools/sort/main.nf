#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSTOOLS_SORT } from '../../../../modules/bustools/sort/main.nf' addParams( options: [:] )

workflow test_bustools_sort {

    input = [ [ id:'test' ], // meta map
            file("https://raw.githubusercontent.com/BUStools/bustools/master/test/output.bus", checkIfExists: true) ]

    BUSTOOLS_SORT ( input )
}
