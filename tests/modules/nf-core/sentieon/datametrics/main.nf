#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_DATAMETRICS } from '../../../../../modules/nf-core/sentieon/datametrics/main.nf'

workflow test_sentieon_datametrics {

    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    fasta = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    SENTIEON_DATAMETRICS ( input, fasta, fai )
}
