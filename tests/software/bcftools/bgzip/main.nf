#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_BGZIP } from '../../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )

workflow test_bcftools_bgzip {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]]

    BCFTOOLS_BGZIP ( input )
}
