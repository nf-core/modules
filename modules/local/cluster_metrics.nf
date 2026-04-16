process CLUSTER_METRICS {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.python_container ?: params.container_py ?: 'snpclustering-py:latest'}"
    publishDir "${params.outdir}/metrics", mode: 'copy'

    input:
    tuple val(meta), path(pca_scores), path(pca_info), path(clusters)
    path cluster_metrics_script

    output:
    tuple val(meta), path("${meta.id}_metrics.tsv"), emit: metrics
    tuple val(meta), path("${meta.id}_k_sweep.csv"), emit: k_sweep
    tuple val(meta), path("${meta.id}_selected.json"), emit: selected
    path "versions.yml", emit: versions

    script:
    def k_min = params.k_min ?: 2
    def k_max = params.k_max ?: 12

    """
    python3 ${cluster_metrics_script} \\
        --features ${pca_scores} \\
        --clusters ${clusters} \\
        --k-min ${k_min} \\
        --k-max ${k_max} \\
        --out-k-sweep ${meta.id}_k_sweep.csv \\
        --out-selected ${meta.id}_selected.json \\
        --out-prefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}