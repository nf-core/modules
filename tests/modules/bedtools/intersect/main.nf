#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_INTERSECT } from '../../../../software/bedtools/intersect/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_intersect {
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
            ]

    BEDTOOLS_INTERSECT ( input )
}
