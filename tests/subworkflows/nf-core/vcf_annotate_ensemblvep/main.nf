#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'

workflow vcf_annotate_ensemblvep {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        []
    ],[
        [ id:'custom_test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
        [
            file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true)
        ]
    ])

    VCF_ANNOTATE_ENSEMBLVEP (
        input,
        [[],[]],
        "WBcel235",
        "caenorhabditis_elegans",
        "108",
        [],
        [],
        5
    )
}

workflow vcf_annotate_ensemblvep_large_chunks {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        []
    ],[
        [ id:'custom_test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true),
        [
            file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test3_vcf_gz_tbi'], checkIfExists: true)
        ]
    ])

    fasta = [
        [id:"fasta"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    VCF_ANNOTATE_ENSEMBLVEP (
        input,
        fasta,
        "WBcel235",
        "caenorhabditis_elegans",
        "108",
        [],
        [],
        100
    )
}
