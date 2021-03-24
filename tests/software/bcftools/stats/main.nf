#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_STATS } from '../../../../software/bcftools/stats/main.nf' addParams( options: [:] )

workflow test_bcftools_stats {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true) ]
            ]

    BCFTOOLS_STATS ( input )
}
