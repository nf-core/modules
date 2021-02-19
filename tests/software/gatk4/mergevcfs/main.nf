#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEVCFS } from '../../../../software/gatk4/mergevcfs/main.nf' addParams( options: [:] )

workflow test_gatk4_mergevcfs {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true),
                file("${launchDir}/tests/data/vcf/test2.vcf.gz", checkIfExists: true)  ] ]

    GATK4_MERGEVCFS ( input )
}
