#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: CLUSTER_VIZ
    Generates PCA, UMAP and t-SNE visualizations colored by cluster
    Author: Donald Baku (athor)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLUSTER_VIZ {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

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
    def prefix = task.ext.prefix ?: out_prefix ?: "${meta.id}"

    """
    python3 ${projectDir}/modules/nf-core/cluster_viz/templates/cluster_viz.py \\
        --features ${features} \\
        --clusters ${clusters} \\
        --pca-scores ${pca_scores} \\
        --out-umap-tsv ${prefix}_umap.tsv \\
        --out-tsne-tsv ${prefix}_tsne.tsv \\
        --out-umap-png ${prefix}_umap.png \\
        --out-tsne-png ${prefix}_tsne.png \\
        --out-pca-png ${prefix}_pca.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
        umap-learn: \$(python3 -c "import umap; print(umap.__version__)" 2>/dev/null || echo 'N/A')
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
