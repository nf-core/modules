#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_FILTER       } from '../../../modules/nf-core/bcftools/filter/main'
include { PLINK2_INDEP_PAIRWISE } from '../../../modules/nf-core/plink2/indeppairwise/main'
include { PLINK2_RECODE_VCF     } from '../../../modules/nf-core/plink2/recodevcf/main'
include { FLASHPCA2             } from '../../../modules/nf-core/flashpca2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SNPCLUSTERING {
    take:
    meta
    vcf
    vcf_index
    maf
    missing

    main:
    versions = Channel.empty()

    BCFTOOLS_FILTER ( vcf.join(vcf_index), maf, missing )
    versions = versions.mix(BCFTOOLS_FILTER.out.versions.first())

    PLINK2_INDEP_PAIRWISE ( BCFTOOLS_FILTER.out.vcf )
    versions = versions.mix(PLINK2_INDEP_PAIRWISE.out.versions.first())

    PLINK2_RECODE_VCF ( PLINK2_INDEP_PAIRWISE.out.pgen )
    versions = versions.mix(PLINK2_RECODE_VCF.out.versions.first())

    FLASHPCA2 ( PLINK2_RECODE_VCF.out.vcf )
    versions = versions.mix(FLASHPCA2.out.versions.first())

    // TODO: qui aggiungeremo KMeans/DBSCAN/plot quando creeremo i moduli local

    emit:
    cluster_labels = Channel.empty()   // placeholder
    metrics        = Channel.empty()   // placeholder
    plots          = Channel.empty()
    versions       = versions
}
