#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GENOMECOV } from '../../../../software/bedtools/genomecov/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_genomecov {
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.bam", checkIfExists: true)
            ]

    BEDTOOLS_GENOMECOV ( input )
}

