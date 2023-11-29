#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WHAMG } from '../../../../modules/nf-core/whamg/main.nf'

workflow test_whamg_bam {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    WHAMG ( input, fasta, fasta_fai )
}
