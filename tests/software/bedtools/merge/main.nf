#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MERGE } from '../../../../software/bedtools/merge/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_merge {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]

    BEDTOOLS_MERGE ( input )
}

