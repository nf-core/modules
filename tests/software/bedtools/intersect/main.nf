#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_INTERSECT } from '../../../../software/bedtools/intersect/main.nf' addParams( options: [:] )

workflow test_bedtools_intersect {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/B.bed", checkIfExists: true) ]

    BEDTOOLS_INTERSECT( input )
}
