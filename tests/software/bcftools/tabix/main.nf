#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_TABIX } from '../../../../software/bcftools/tabix/main.nf' addParams( options: [:] )

workflow test_bcftools_tabix {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/vcf/test.vcf.gz", checkIfExists: true) ] ]

    BCFTOOLS_TABIX ( input )
}
