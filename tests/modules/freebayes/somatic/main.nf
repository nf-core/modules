#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREEBAYES_SOMATIC } from '../../../../modules/freebayes/somatic/main.nf' addParams( options: [:] )

workflow test_freebayes {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai     = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    targets = []
    samples = []

    FREEBAYES_SOMATIC ( input, fasta, fai, targets, samples)
}

workflow test_freebayes_intervals {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram_crai'], checkIfExists: true)
            ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai     = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    targets = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    samples = []

    FREEBAYES_SOMATIC ( input, fasta, fai, targets, samples)
}
