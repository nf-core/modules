#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GENOMECOV } from '../../../../software/bedtools/genomecov/main.nf' addParams( options: [:] )

workflow test_bedtools_genomecov {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bam/sarscov2_paired_aln.bam", checkIfExists: true) ]

    BEDTOOLS_GENOMECOV ( input )
}

