#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABRA2_REALIGN as ABRA2_REALIGN_SINGLE_END } from '../../../../../modules/nf-core/abra2/realign/main.nf'
include { ABRA2_REALIGN as ABRA2_REALIGN_ONE } from '../../../../../modules/nf-core/abra2/realign/main.nf'
include { ABRA2_REALIGN as ABRA2_REALIGN_PAIR} from '../../../../../modules/nf-core/abra2/realign/main.nf'
include { ABRA2_REALIGN as ABRA2_REALIGN_TRIO} from '../../../../../modules/nf-core/abra2/realign/main.nf'

workflow test_abra2_realign_single_end_sample {

    input = [
        [ id:'test-alone', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    genome = [
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    targets = file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)

    input_empty = [[:],[],[]]

    ABRA2_REALIGN_SINGLE_END ( input, input_empty, input_empty, genome, targets )
}

workflow test_abra2_realign_one_sample {

    input = [
        [ id:'test-alone', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    genome = [
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    targets = file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)

    input_empty = [[:],[],[]]

    ABRA2_REALIGN_ONE ( input, input_empty, input_empty, genome, targets )
}

workflow test_abra2_realign_two_samples {

    input1 = [
        [ id:'test-paired-samp', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    input2 = [
        [ id:'test2-paired-samp', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    input_empty = [[:],[],[]]

    genome = [
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    targets = file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)

    ABRA2_REALIGN_PAIR ( input1, input2,     input_empty, genome, targets )
}

workflow test_abra2_realign_three_samples {

    input1 = [
        [ id:'test-paired-samp', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    input2 = [
        [ id:'test2-paired-samp', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    input_rna = [
        [ id:'test_rna', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_rna_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_rna_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    input_empty = [[:],[],[]]

    genome = [
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    targets = file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)

    ABRA2_REALIGN_TRIO ( input1, input2,     input_rna,   genome, targets )
}
