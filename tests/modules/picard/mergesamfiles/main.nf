#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MERGESAMFILES } from '../../../../modules/picard/mergesamfiles/main.nf'

workflow test_picard_mergesamfiles {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true), ]
            ]

    PICARD_MERGESAMFILES ( input )
}
