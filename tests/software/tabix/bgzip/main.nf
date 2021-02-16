#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_BGZIP as BGZIP_VCF } from '../../../../software/tabix/bgzip/main.nf' addParams( options: [:] )

workflow test_tabix_bgzip_vcf {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]]

    BGZIP_VCF ( input )
}
