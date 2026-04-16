#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PLINK2_VCF }         from '../../../modules/nf-core/plink2/vcf/main'
include { PLINK2_PGEN2BED }    from '../../../modules/local/plink2_pgen2bed'
include { PCA_FLASHPCA }       from '../../../modules/local/pca'
include { CLUSTERING }         from '../../../modules/local/clustering'
include { CLUSTER_METRICS }    from '../../../modules/local/cluster_metrics'
include { CLUSTER_VIZ }        from '../../../modules/local/cluster_viz'

workflow SNPCLUSTERING {

    take:
    vcf_ch

    main:
    def versions_ch = Channel.empty()

    PLINK2_VCF(vcf_ch)
    versions_ch = versions_ch.mix(PLINK2_VCF.out.versions.ifEmpty([]))

    ch_pgen_triple = PLINK2_VCF.out.pgen
        .join(PLINK2_VCF.out.pvar)
        .join(PLINK2_VCF.out.psam)

    PLINK2_PGEN2BED(ch_pgen_triple)
    versions_ch = versions_ch.mix(PLINK2_PGEN2BED.out.versions.ifEmpty([]))

    ch_bed_triple = PLINK2_PGEN2BED.out.bed
        .join(PLINK2_PGEN2BED.out.bim)
        .join(PLINK2_PGEN2BED.out.fam)

    parser_script_ch = Channel.value(
        file("${projectDir}/subworkflows/nf-core/snpclustering/scripts/flashpca_outpc_to_tsv.py", checkIfExists: true)
    )

    PCA_FLASHPCA(ch_bed_triple, parser_script_ch)
    ch_pca = PCA_FLASHPCA.out.pca
    versions_ch = versions_ch.mix(PCA_FLASHPCA.out.versions.ifEmpty([]))

    ch_pca_for_clustering = ch_pca.map { meta, features, scaled, pca_scores, pca_info ->
        tuple(meta, pca_scores, pca_info)
    }
    
    clustering_script_ch = Channel.value(
    file("${projectDir}/subworkflows/nf-core/snpclustering/scripts/clustering.py", checkIfExists: true)
)
    CLUSTERING(ch_pca_for_clustering, clustering_script_ch)

    ch_for_metrics_viz = ch_pca_for_clustering.join(CLUSTERING.out.clusters)
    cluster_metrics_script_ch = Channel.value(
    file("${projectDir}/subworkflows/nf-core/snpclustering/scripts/cluster_metrics.py", checkIfExists: true)
)
    CLUSTER_METRICS(ch_for_metrics_viz, cluster_metrics_script_ch)

    cluster_viz_script_ch = Channel.value(
    file("${projectDir}/subworkflows/nf-core/snpclustering/scripts/cluster_viz.py", checkIfExists: true)
)
    CLUSTER_VIZ(ch_for_metrics_viz, cluster_viz_script_ch)

    versions_ch = versions_ch
        .mix(CLUSTER_METRICS.out.versions.ifEmpty([]))
        .mix(CLUSTER_VIZ.out.versions.ifEmpty([]))
        .unique()

    emit:
    plink_bed      = ch_bed_triple
    pca            = ch_pca
    cluster_labels = CLUSTERING.out.clusters
    metrics        = CLUSTER_METRICS.out.metrics
    plots          = CLUSTER_VIZ.out.plots
    versions       = versions_ch
}