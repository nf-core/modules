#!/usr/bin/env nextflow



include { BEDTOOLS_MERGE } from '../../../../modules/bedtools/merge/main.nf'

workflow test_bedtools_merge {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]

    BEDTOOLS_MERGE ( input )
}

