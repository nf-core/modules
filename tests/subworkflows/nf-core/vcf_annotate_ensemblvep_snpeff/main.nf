#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_DOWNLOAD            } from '../../../../modules/nf-core/ensemblvep/download/main'
include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main'

vep_cache_version = "110"
vep_genome = "WBcel235"
vep_species = "caenorhabditis_elegans"
vep_cache_input = Channel.of([[id:"${vep_cache_version}_${vep_genome}"], vep_genome, vep_species, vep_cache_version])

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

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        [[],[]],
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
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

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        [[],[]],
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        [],
        "WBcel235.99",
        [],
        ["ensemblvep", "snpeff"],
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

    fasta = Channel.value([
        [id:"fasta"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ])

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        fasta,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
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

    fasta = Channel.value([
        [id:"fasta"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ])

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        fasta,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        [],
        [],
        [],
        ["ensemblvep"],
        []
    )
}
