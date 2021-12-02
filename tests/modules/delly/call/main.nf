#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DELLY_CALL } from '../../../../modules/delly/call/main.nf'

workflow test_delly_call {
    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
                ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DELLY_CALL ( input, fasta, fai )
}
