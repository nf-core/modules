#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { JASMINESV } from '../../../../modules/nf-core/jasminesv/main.nf'

workflow test_jasminesv_minimum {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf'], checkIfExists: true)
        ],
        [],
        []
    ]

    fasta = []
    fasta_fai = []
    chr_norm = []

    JASMINESV ( input, fasta, fasta_fai, chr_norm )
}

workflow test_jasminesv_bgzipped_input {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true)
        ],
        [],
        []
    ]

    fasta = []
    fasta_fai = []
    chr_norm = []

    JASMINESV ( input, fasta, fasta_fai, chr_norm )
}

workflow test_jasminesv_iris {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf'], checkIfExists: true)
        ],
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
        ],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    chr_norm = []

    JASMINESV ( input, fasta, fasta_fai, chr_norm )
}

workflow test_jasminesv_all_inputs {

    input = Channel.of(
        [
            [ id:'test', single_end:false ], // meta map
            [
                file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf'], checkIfExists: true)
            ],
            [
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
            ],
            []
        ]
    )

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    chr_norm = Channel.of("chr21 21", "chr22 22").collectFile(name:"chr_norm.txt", newLine:true)

    JASMINESV ( input, fasta, fasta_fai, chr_norm )
}
