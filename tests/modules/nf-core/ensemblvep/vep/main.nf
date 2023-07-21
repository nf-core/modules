#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_DEFAULT   } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_JSON      } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_TAB       } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF       } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF_BGZIP } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF_GZIP  } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_CUSTOM    } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_DOWNLOAD                        } from '../../../../../modules/nf-core/ensemblvep/download/main.nf'

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

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_JSON ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, fasta, [] )
}

workflow test_ensemblvep_vep_fasta_tab {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    fasta = [
        [ id: 'fasta' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    ENSEMBLVEP_VEP_TAB ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, fasta, [] )
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

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, fasta, [] )
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

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF_BGZIP ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, fasta, [] )
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

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF_GZIP ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, fasta, [] )
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

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_DEFAULT ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, fasta, [] )
}

workflow test_ensemblvep_vep_no_fasta {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ]

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_DEFAULT ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, [[], []], [] )
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

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF_BGZIP ( input, "WBcel235", "caenorhabditis_elegans", "110", cache, fasta, [] )
}
