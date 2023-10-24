#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARDICTJAVA } from '../../../../modules/nf-core/vardictjava/main.nf'

workflow test_vardictjava_bam {
    
    bam_input_ch = Channel.value([
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ])

    fasta = [
        [id:"ref"],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [id:"ref"],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VARDICTJAVA(
        bam_input_ch,
        fasta,
        fasta_fai
    )
}

workflow test_vardictjava_paired_bam {
    
    bam_input_ch = Channel.value([
        [ id:'test' ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
        ],
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
        ],
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ])

    fasta = [
        [id:"ref"],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [id:"ref"],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VARDICTJAVA(
        bam_input_ch,
        fasta,
        fasta_fai
    )
}
