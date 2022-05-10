#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS } from '../../../../modules/gatk4/collectreadcounts/main.nf'

workflow test_gatk4_collectreadcounts {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_COLLECTREADCOUNTS ( input, fasta, fai, dict )
}
