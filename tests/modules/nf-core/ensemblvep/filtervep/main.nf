#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_DOWNLOAD  } from '../../../../../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_FILTERVEP as ENSEMBLVEP_FILTERVEP_TAB } from '../../../../../modules/nf-core/ensemblvep/filtervep/main.nf'
include { ENSEMBLVEP_FILTERVEP as ENSEMBLVEP_FILTERVEP_VCF } from '../../../../../modules/nf-core/ensemblvep/filtervep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_TAB } from '../../../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF } from '../../../../../modules/nf-core/ensemblvep/vep/main'

vep_cache_version = "110"
vep_genome = "WBcel235"
vep_species = "caenorhabditis_elegans"
vep_cache_input = Channel.of([[id:"${vep_cache_version}_${vep_genome}"], vep_genome, vep_species, vep_cache_version])

workflow test_ensemblvep_filtervep_tab {
    
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ])

    fasta = Channel.value([
        [id:"fasta"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ])

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_TAB ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )

    ENSEMBLVEP_FILTERVEP_TAB ( ENSEMBLVEP_VEP_TAB.out.tab, [] )
}

workflow test_ensemblvep_filtervep_vcf {
    
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        []
    ])

    fasta = Channel.value([
        [id:"fasta"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ])

    ENSEMBLVEP_DOWNLOAD(vep_cache_input)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    ENSEMBLVEP_VEP_VCF ( input, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, [] )

    ENSEMBLVEP_FILTERVEP_VCF ( ENSEMBLVEP_VEP_VCF.out.vcf, [] )
}
