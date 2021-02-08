#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_CONSENSUS } from '../../../../software/bcftools/consensus/main.nf' addParams( options: [:] )

workflow test_bcftools_consensus {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/vcf/test.vcf.gz", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/vcf/test.vcf.gz.tbi", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/vcf/test.consensus.fa", checkIfExists: true) ] ]

    BCFTOOLS_CONSENSUS ( input )
}
