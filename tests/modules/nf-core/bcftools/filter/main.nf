#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_FILTER } from '../../../../../modules/nf-core/bcftools/filter/main.nf'

workflow test_bcftools_filter_vcf_gz {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    BCFTOOLS_FILTER ( input )
}

workflow test_bcftools_filter_bcf {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    BCFTOOLS_FILTER ( input )
}

