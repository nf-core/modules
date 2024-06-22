#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_DOWNLOAD            } from '../../../../modules/nf-core/ensemblvep/download/main'
include { SNPEFF_DOWNLOAD                } from '../../../../modules/nf-core/snpeff/download/main'
include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main'

snpeff_cache_version = "105"
snpeff_genome = "WBcel235"
snpeff_cache_input = Channel.of([[id:"${snpeff_genome}.${snpeff_cache_version}"], snpeff_genome, snpeff_cache_version])
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

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }.first()

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

    SNPEFF_DOWNLOAD(snpeff_cache_input)

    snpeff_cache = SNPEFF_DOWNLOAD.out.cache.first()

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        [[],[]],
        [],
        [],
        [],
        [],
        [],
        "${snpeff_genome}.${snpeff_cache_version}",
        snpeff_cache,
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
    SNPEFF_DOWNLOAD(snpeff_cache_input)

    snpeff_cache = SNPEFF_DOWNLOAD.out.cache.first()
    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }.first()

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF (
        input,
        [[],[]],
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        [],
        "${snpeff_genome}.${snpeff_cache_version}",
        snpeff_cache,
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

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }.first()

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

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }.first()

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
