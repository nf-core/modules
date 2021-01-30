#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MERGE } from '../../../../software/bedtools/merge/main.nf' addParams( options: [:] )


workflow test_bedtools_merge {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]

    BEDTOOLS_MERGE( input )
}


