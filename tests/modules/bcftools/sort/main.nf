#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_SORT } from '../../../../modules/bcftools/sort/main.nf'

workflow test_bcftools_sort {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    BCFTOOLS_SORT ( input )
}
