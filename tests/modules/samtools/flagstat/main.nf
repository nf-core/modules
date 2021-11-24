#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT } from '../../../../modules/samtools/flagstat/main.nf'

workflow test_samtools_flagstat {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    SAMTOOLS_FLAGSTAT ( input )
}
