#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MANTA_GERMLINE } from '../../../../modules/manta/germline/main.nf' addParams( options: [:] )

workflow test_manta_germline {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true)
    ]

    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai     = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    bed     = []
    bed_tbi = []

    MANTA_GERMLINE ( input, fasta, fai, bed, bed_tbi )
}

workflow test_manta_germline_target_bed {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_crai'], checkIfExists: true)
    ]

    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai     = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    bed     = file(params.test_data['homo_sapiens']['genome']['genome_bed_gz'], checkIfExists: true)
    bed_tbi = file(params.test_data['homo_sapiens']['genome']['genome_bed_gz_tbi'], checkIfExists: true)

    MANTA_GERMLINE ( input, fasta, fai, bed, bed_tbi )
}
