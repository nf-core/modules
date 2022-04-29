#!/usr/bin/env nextflow



//keep --no-verson argument, otherwise md5 will change on each execution
include { BCFTOOLS_FILTER } from '../../../../modules/bcftools/filter/main.nf'

workflow test_bcftools_filter {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]

    BCFTOOLS_FILTER ( input )
}
