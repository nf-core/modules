#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main.nf'

workflow vcf_annotate_ensemblvep_snpeff_vep {
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

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        [[],[]],
        "WBcel235",
        "caenorhabditis_elegans",
        "108",
        [],
        [],
        [],
        [],
        ["ensemblvep"],
        5
    )
}

workflow vcf_annotate_ensemblvep_snpeff_snpeff {
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

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        [[],[]],
        [],
        [],
        [],
        [],
        [],
        "WBcel235.99",
        [],
        ["snpeff"],
        5
    )
}

workflow vcf_annotate_ensemblvep_snpeff_both {
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

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        [[],[]],
        "WBcel235",
        "caenorhabditis_elegans",
        "108",
        [],
        [],
        "WBcel235.99",
        [],
        ["snpeff", "ensemblvep"],
        5
    )
}

workflow vcf_annotate_ensemblvep_snpeff_large_chunks {
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

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        fasta,
        "WBcel235",
        "caenorhabditis_elegans",
        "108",
        [],
        [],
        [],
        [],
        ["ensemblvep"],
        100
    )
}

workflow vcf_annotate_ensemblvep_snpeff_no_scatter {
    input = Channel.of([
        [ id:'test1' ], // meta map
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

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        fasta,
        "WBcel235",
        "caenorhabditis_elegans",
        "108",
        [],
        [],
        [],
        [],
        ["ensemblvep"],
        []
    )
}
