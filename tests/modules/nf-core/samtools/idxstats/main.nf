#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_IDXSTATS } from '../../../../modules/samtools/idxstats/main.nf'

workflow test_samtools_idxstats {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    SAMTOOLS_IDXSTATS ( input )
}
