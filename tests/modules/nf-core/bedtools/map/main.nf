#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MAP                           } from '../../../../../modules/nf-core/bedtools/map/main.nf'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_BAM } from '../../../../../modules/nf-core/bedtools/map/main.nf'

workflow test_bedtools_map {
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
    ]


    BEDTOOLS_MAP ( input, [[:],[]] )
}

workflow test_bedtools_map_bam {
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    BEDTOOLS_MAP_BAM ( input, [[:],[]] )
}
