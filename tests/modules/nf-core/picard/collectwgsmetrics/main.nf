#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTWGSMETRICS } from '../../../../../modules/nf-core/picard/collectwgsmetrics/main.nf'

workflow test_picard_collectwgsmetrics {
    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    PICARD_COLLECTWGSMETRICS ( input, fasta )
}
