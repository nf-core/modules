#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTHSMETRICS } from '../../../../modules/picard/collecthsmetrics/main.nf'

workflow test_picard_collecthsmetrics {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    bait_intervals = file(params.test_data['sarscov2']['genome']['baits_interval_list'], checkIfExists: true)
    target_intervals = file(params.test_data['sarscov2']['genome']['targets_interval_list'], checkIfExists: true)

    PICARD_COLLECTHSMETRICS ( input, fasta, fai, bait_intervals, target_intervals )
}
