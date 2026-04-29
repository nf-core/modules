#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: CLUSTERING
    Unified PLINK2 PCA + Clustering (KMeans / DBSCAN)
    Author: Donald Baku (athor)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLUSTERING {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(eigenvec)
    val algorithm
    val n_clusters
    val dbscan_eps
    val dbscan_min_samples

    output:
    tuple val(meta), path("*_clusters.csv")         , emit: clusters
    tuple val(meta), path("*_clustering_info.json") , emit: info, optional: true
    path "versions.yml"                             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 ${projectDir}/modules/nf-core/clustering/templates/plink2_clustering.py \\
        --eigenvec ${eigenvec} \\
        --algorithm ${algorithm} \\
        --k ${n_clusters} \\
        --dbscan-eps ${dbscan_eps} \\
        --dbscan-min-samples ${dbscan_min_samples} \\
        --out-prefix ${prefix} \\
        --id-mode iid

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
