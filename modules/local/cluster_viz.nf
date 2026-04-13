process CLUSTER_VIZ {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.python_container ?: params.container_py ?: 'snpclustering-py:latest'}"
    publishDir "${params.outdir}/viz", mode: 'copy'

    input:
    tuple val(meta), path(pca_scores), path(pca_info), path(clusters)
    path cluster_viz_script

    output:
    tuple val(meta), path("${meta.id}_umap.png"),
                     path("${meta.id}_tsne.png"),
                     path("${meta.id}_pca.png"), emit: plots
    tuple val(meta), path("${meta.id}_umap.tsv"), emit: umap_tsv
    tuple val(meta), path("${meta.id}_tsne.tsv"), emit: tsne_tsv
    path "versions.yml", emit: versions

    script:
    def umap_n = params.viz_umap_neighbors ?: 15
    def umap_d = params.viz_umap_min_dist ?: 0.1
    def tsne_p = params.viz_perplexity ?: 30
    def tsne_i = params.viz_tsne_iter ?: 1000

    """
    export NUMBA_DISABLE_JIT=1

    python3 ${cluster_viz_script} \\
        --features ${pca_scores} \\
        --clusters ${clusters} \\
        --pca-scores ${pca_scores} \\
        --umap-neighbors ${umap_n} \\
        --umap-min-dist ${umap_d} \\
        --tsne-perplexity ${tsne_p} \\
        --tsne-iter ${tsne_i} \\
        --out-umap-tsv ${meta.id}_umap.tsv \\
        --out-tsne-tsv ${meta.id}_tsne.tsv \\
        --out-umap-png ${meta.id}_umap.png \\
        --out-tsne-png ${meta.id}_tsne.png \\
        --out-pca-png ${meta.id}_pca.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}