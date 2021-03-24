#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_FILTER } from '../../../../software/bcftools/filter/main.nf' addParams( options: ['args': '--no-version'] )

workflow test_bcftools_filter {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true) ]
            ]

    BCFTOOLS_FILTER ( input )
}
