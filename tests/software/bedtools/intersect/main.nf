#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_INTERSECT } from '../../../../software/bedtools/intersect/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_intersect {
    input = [ [ id:'test' ],
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test2.bed", checkIfExists: true)
            ]

    BEDTOOLS_INTERSECT ( input )
}
