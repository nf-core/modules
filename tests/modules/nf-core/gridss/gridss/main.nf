#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GRIDSS_GRIDSS } from '../../../../../modules/nf-core/gridss/gridss/main.nf'
include { BWA_INDEX     } from '../../../../../modules/nf-core/bwa/index/main.nf'

workflow test_gridss_gridss_bwa {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        []
    ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    bwa_index = BWA_INDEX(fasta).index

    GRIDSS_GRIDSS(
        input,
        fasta,
        fasta_fai,
        bwa_index
    )
}

workflow test_gridss_gridss_no_bwa {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        []
    ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    GRIDSS_GRIDSS(
        input,
        fasta,
        fasta_fai,
        [[], []]
    )
}

workflow test_gridss_gridss_assembly {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
    ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    GRIDSS_GRIDSS(
        input,
        fasta,
        fasta_fai,
        [[], []]
    )
}

workflow test_gridss_gridss_cram {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        []
    ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    GRIDSS_GRIDSS(
        input,
        fasta,
        fasta_fai,
        [[], []]
    )
}
