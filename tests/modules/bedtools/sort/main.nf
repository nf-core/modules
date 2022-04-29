#!/usr/bin/env nextflow



include { BEDTOOLS_SORT } from '../../../../modules/bedtools/sort/main.nf'

workflow test_bedtools_sort {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]

    BEDTOOLS_SORT ( input, "testext" )
}
