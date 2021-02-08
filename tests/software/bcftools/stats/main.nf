#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_STATS } from '../../../../software/bcftools/stats/main.nf' addParams( options: [:] )

workflow test_bcftools_stats {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]]

    BCFTOOLS_STATS ( input )
}
