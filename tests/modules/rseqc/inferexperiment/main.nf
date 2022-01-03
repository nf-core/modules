#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEQC_INFEREXPERIMENT }   from '../../../../modules/rseqc/inferexperiment/main.nf'

workflow test_rseqc_inferexperiment {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    RSEQC_INFEREXPERIMENT ( input, bed )
}
