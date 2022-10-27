#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_FILTER } from "$moduleDir/modules/nf-core/bcftools/filter/main.nf"

workflow test_bcftools_filter {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    BCFTOOLS_FILTER ( input )
}
