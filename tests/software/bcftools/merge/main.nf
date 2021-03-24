#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_MERGE } from '../../../../software/bcftools/merge/main.nf' addParams( options: ['args': '--force-samples --no-version'] )

workflow test_bcftools_merge {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test2.vcf.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test3.vcf.gz", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test2.vcf.gz.tbi", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test3.vcf.gz.tbi", checkIfExists: true) ]
            ]

    BCFTOOLS_MERGE ( input )
}
