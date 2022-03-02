#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEQC_TIN } from '../../../../modules/rseqc/tin/main.nf'

workflow test_rseqc_tin {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    RSEQC_TIN ( input, bed )
}
