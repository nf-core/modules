#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: CLUSTER_VIZ
    Generates PCA, UMAP and t-SNE visualizations colored by cluster
    Author: Donald Baku (author)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLUSTER_VIZ {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/clustervisualiation:dev' :
    'quay.io/nf-core/clustervisualiation:dev' }"
    input:
    tuple val(meta), path(features), path(clusters), path(pca_scores)
    val out_prefix

    output:
    tuple val(meta), path("*_umap.tsv")   , emit: umap
    tuple val(meta), path("*_tsne.tsv")   , emit: tsne
    tuple val(meta), path("*_umap.png")   , emit: umap_png
    tuple val(meta), path("*_tsne.png")   , emit: tsne_png
    tuple val(meta), path("*_pca.png")    , emit: pca_png
    path "versions.yml"                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'cluster_viz.py'

}
