#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_INTERSECT } from '../../../../modules/bedtools/intersect/main.nf'

workflow test_bedtools_intersect {
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
    ]

    extension = 'bed'

    BEDTOOLS_INTERSECT ( input, extension )
}

workflow test_bedtools_intersect_bam {
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true)
    ]

    extension = 'bam'

    BEDTOOLS_INTERSECT ( input, extension )
}
