#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SNPCLUSTERING
    Author: Donald
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PLINK2_VCF } from '../../../modules/nf-core/plink2/vcf/main'
include { PLINK2_PCA } from '../../../modules/nf-core/plink2/pca/main'
include { CLUSTERING } from '../../../modules/nf-core/clustering/main'
include { CLUSTER_METRICS } from '../../../modules/nf-core/cluster_metrics/main'
include { CLUSTER_VIZ } from '../../../modules/nf-core/cluster_viz/main'

workflow SNPCLUSTERING {
    take:
    ch_vcf // channel: [ val(meta), path(vcf) ]

    main:
    ch_versions = Channel.empty()

    // 1. VCF → PGEN
    PLINK2_VCF(ch_vcf)
    ch_versions = ch_versions.mix(PLINK2_VCF.out.versions)

    ch_pgen = PLINK2_VCF.out.pgen
        .join(PLINK2_VCF.out.pvar)
        .join(PLINK2_VCF.out.psam)

    // 2. PCA
    PLINK2_PCA(
        ch_pgen.map { meta, pgen, pvar, psam ->
            [ meta, params.npcs ?: 10, params.use_approx ?: true, pgen, psam, pvar ]
        }
    )
    ch_versions = ch_versions.mix(PLINK2_PCA.out.versions)

    // 3. Clustering
    CLUSTERING(
        PLINK2_PCA.out.evecfile,
        params.algorithm ?: 'kmeans',
        params.n_clusters ?: 3,
        params.dbscan_eps ?: 0.5,
        params.dbscan_min_samples ?: 5
    )
    ch_versions = ch_versions.mix(CLUSTERING.out.versions)

    // 4. Join per downstream
    ch_for_downstream = PLINK2_PCA.out.evecfile
        .join(CLUSTERING.out.clusters, by: 0)

    // 5. Metrics + Viz (se vuoi includerli)
    // CLUSTER_METRICS(ch_for_downstream)
    // ch_versions = ch_versions.mix(CLUSTER_METRICS.out.versions)

    // CLUSTER_VIZ(ch_for_downstream)
    // ch_versions = ch_versions.mix(CLUSTER_VIZ.out.versions)

    emit:
    pca_eigenvec = PLINK2_PCA.out.evecfile
    pca_eigenval = PLINK2_PCA.out.evfile
    clusters     = CLUSTERING.out.clusters
    // metrics      = CLUSTER_METRICS.out.metrics
    // plots        = CLUSTER_VIZ.out.plots
    versions     = ch_versions
}
