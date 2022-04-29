#!/usr/bin/env nextflow



include { BEDTOOLS_MAKEWINDOWS } from '../../../../modules/bedtools/makewindows/main.nf'

workflow test_bedtools_makewindows {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    BEDTOOLS_MAKEWINDOWS ( input, true )
}
