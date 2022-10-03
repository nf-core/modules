#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fragment_size = 1000

include { DEEPTOOLS_PLOTFINGERPRINT } from '../../../../modules/deeptools/plotfingerprint/main.nf'

workflow test_deeptools_plotfingerprint {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    DEEPTOOLS_PLOTFINGERPRINT ( input )
}
