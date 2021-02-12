#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_MERGE } from '../../../../software/bcftools/merge/main.nf' addParams( options: ['args': '--force-samples --no-version'] )

workflow test_bcftools_merge {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/vcf/test_a.merge.vcf.gz", checkIfExists: true),
                file("${launchDir}/tests/data/vcf/test_b.merge.vcf.gz", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/vcf/test_a.merge.vcf.gz.tbi", checkIfExists: true),
                file("${launchDir}/tests/data/vcf/test_b.merge.vcf.gz.tbi", checkIfExists: true) ]]

    BCFTOOLS_MERGE ( input )
}
