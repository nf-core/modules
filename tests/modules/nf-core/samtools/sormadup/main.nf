#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_SORMADUP } from '../../../../../modules/nf-core/samtools/sormadup/main.nf'

workflow test_samtools_sormadup {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        ]
    ]

    SAMTOOLS_SORMADUP ( input, [[:],[]] )
}
