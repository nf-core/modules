#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_PRINTREADS } from '../../../../../modules/nf-core/gatk4/printreads/main.nf'

workflow test_gatk4_printreads_bam {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)

    ]

    fasta = [
        [], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]

    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_PRINTREADS ( input, fasta, fai, dict )
}

workflow test_gatk4_printreads_cram {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)

    ]

    fasta = [
        [], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
        ]

    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_PRINTREADS ( input, fasta, fai, dict )
}
