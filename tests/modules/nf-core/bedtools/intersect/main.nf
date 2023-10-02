#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_INTERSECT                           } from '../../../../../modules/nf-core/bedtools/intersect/main.nf'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_BAM } from '../../../../../modules/nf-core/bedtools/intersect/main.nf'

workflow test_bedtools_intersect {
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
    ]


    BEDTOOLS_INTERSECT ( input, [[:],[]] )
}

workflow test_bedtools_intersect_bam {
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true)
    ]

    BEDTOOLS_INTERSECT_BAM ( input, [[:],[]] )
}
