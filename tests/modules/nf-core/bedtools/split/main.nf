#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SPLIT } from '../../../../modules/bedtools/split/main.nf'

workflow test_bedtools_split {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ]

    number_of_files = 2

    BEDTOOLS_SPLIT ( input, number_of_files )
}
