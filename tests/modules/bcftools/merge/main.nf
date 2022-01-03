#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_MERGE } from '../../../../modules/bcftools/merge/main.nf'

workflow test_bcftools_merge {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true) ]
            ]

    BCFTOOLS_MERGE ( input )
}
