#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_DEFAULT   } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_JSON      } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_TAB       } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF       } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF_BGZIP } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF_GZIP  } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_CUSTOM    } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_DOWNLOAD                        } from '../../../../../modules/nf-core/ensemblvep/download/main'

vep_cache_version = "110"
vep_genome = "WBcel235"
vep_species = "caenorhabditis_elegans"
vep_cache_input = Channel.of([[id:"${vep_cache_version}_${vep_genome}"], vep_genome, vep_species, vep_cache_version])

workflow test_ensemblvep_vep_fasta_json {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_JSON ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )
}

workflow test_ensemblvep_vep_fasta_tab {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_VEP_TAB ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )
}

workflow test_ensemblvep_vep_fasta_vcf {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )
}

workflow test_ensemblvep_vep_fasta_vcf_bgzip {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF_BGZIP ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )
}

workflow test_ensemblvep_vep_fasta_vcf_gzip {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF_GZIP ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )
}

workflow test_ensemblvep_vep_fasta {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_DEFAULT ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )
}

workflow test_ensemblvep_vep_no_fasta {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_DEFAULT ( input, vep_genome, vep_species, vep_cache_version, vep_cache, [[], []], [] )
}

workflow test_ensemblvep_vep_fasta_custom {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
            [ file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test3_vcf'], checkIfExists: true)]
    ]

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF_BGZIP ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )
}
