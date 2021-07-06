#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SORT } from '../../../../software/bedtools/sort/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_sort {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]

    BEDTOOLS_SORT ( input )
}
