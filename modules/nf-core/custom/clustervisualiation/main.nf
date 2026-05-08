process CUSTOM_CLUSTERVISUALIATION {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c00b83d40a02e4ed2833ebf0d38635602231a21764eff0d30ed16885e5c02445/data' :
    'community.wave.seqera.io/library/matplotlib_pandas_python_scikit-learn_umap-learn:2c4aaf377be5cd4a' }"

    input:
    tuple val(meta), path(features), path(clusters), path(pca_scores)

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
