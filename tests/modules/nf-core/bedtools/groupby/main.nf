#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GROUPBY } from '../../../../../modules/nf-core/bedtools/groupby/main.nf'

workflow test_bedtools_groupby {

    input = [
        [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
    ]
    BEDTOOLS_GROUPBY ( input, 5 )
}
