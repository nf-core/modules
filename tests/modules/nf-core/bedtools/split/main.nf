#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SPLIT } from '../../../../../modules/nf-core/bedtools/split/main.nf'

workflow test_bedtools_split {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true),
        2
    ]

    BEDTOOLS_SPLIT ( input )
}
