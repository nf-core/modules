#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEQC_JUNCTIONSATURATION }   from '../../../../modules/rseqc/junctionsaturation/main.nf'

workflow test_rseqc_junctionsaturation {
    input = [
        [ id:'test', single_end: false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    RSEQC_JUNCTIONSATURATION ( input, bed )
}
