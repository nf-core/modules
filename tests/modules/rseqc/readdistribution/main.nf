#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEQC_READDISTRIBUTION }   from '../../../../modules/rseqc/readdistribution/main.nf'

workflow test_rseqc_readdistribution {
    input = [ [ id:'test', single_end: false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    bed = file(params.test_data['sarscov2']['genome']['test_bed12'], checkIfExists: true)

    RSEQC_READDISTRIBUTION ( input, bed )
}
