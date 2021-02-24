#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SORT } from '../../../../software/bedtools/sort/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_sort {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true) ]

    BEDTOOLS_SORT ( input )
}
