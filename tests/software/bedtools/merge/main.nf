#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MERGE } from '../../../../software/bedtools/merge/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_merge {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true) ]

    BEDTOOLS_MERGE ( input )
}

