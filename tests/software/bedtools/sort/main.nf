#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SORT } from '../../../../software/bedtools/sort/main.nf' addParams( options: [:] )

workflow test_bedtools_sort {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]

    BEDTOOLS_SORT ( input )
}
