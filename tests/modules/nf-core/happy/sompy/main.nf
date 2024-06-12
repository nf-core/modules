#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPPY_SOMPY } from '../../../../../modules/nf-core/happy/sompy/main.nf'

workflow test_sompy {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        [],
        []
    ]

    fasta = [
        [ id:'test1' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [ id:'test1' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    HAPPY_SOMPY (
        input,
        fasta,
        fasta_fai,
        [[],[]],
        [[],[]],
        [[],[]]
    )
}
workflow test_sompy_bam {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        [],
        []
    ]

    fasta = [
        [ id:'test1' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta_fai = [
        [ id:'test1' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    bam = [
        [id:'test2'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists:true)
    ]

    HAPPY_SOMPY (
        input,
        fasta,
        fasta_fai,
        [[],[]],
        [[],[]],
        bam
    )
}
