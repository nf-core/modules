#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_BGZIP } from '../../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX } from '../../../../software/bcftools/tabix/main.nf' addParams( options: [:] )

workflow test_bcftools_tabix {

    def input = []
    input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    BCFTOOLS_BGZIP ( input )
    BCFTOOLS_TABIX ( BCFTOOLS_BGZIP.out.vcf )
}
