#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIP } from '../../../../software/tabix/bgzip/main.nf' addParams( options: [:] )

workflow test_tabix_bgzip {
    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf", checkIfExists: true) ]]

    TABIX_BGZIP ( input )
}
