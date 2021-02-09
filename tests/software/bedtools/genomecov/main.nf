#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GENOMECOV } from '../../../../software/bedtools/genomecov/main.nf' addParams( options: [:] )

workflow test_bedtools_genomecov {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bam/test.paired_end.name.sorted.bam", checkIfExists: true) ]

    BEDTOOLS_GENOMECOV ( input )
}

