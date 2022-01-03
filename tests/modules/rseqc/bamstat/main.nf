#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEQC_BAMSTAT }   from '../../../../modules/rseqc/bamstat/main.nf'

workflow test_rseqc_bamstat {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    RSEQC_BAMSTAT ( input )
}
