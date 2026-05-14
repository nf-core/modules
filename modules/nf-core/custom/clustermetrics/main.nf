process CUSTOM_CLUSTERMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/25129a5258522a434c386b800d3e2e3e6dc72d8a1171b7b10f21df3488526795/data' :
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:169e228afc7d3686' }"

    input:
    tuple val(meta), path(features), path(clusters)

    output:
    tuple val(meta), path("*.metrics.tsv")  , emit: metrics
    tuple val(meta), path("*.k_sweep.csv")  , emit: k_sweep
    tuple val(meta), path("*.selected.json"), emit: selected
    tuple val(meta), path("*.png")          , emit: plots, optional: true
    path "versions.yml"                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'cluster_metrics.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.metrics.tsv
    touch ${prefix}.k_sweep.csv
    touch ${prefix}.selected.json
    touch ${prefix}.elbow.png
    touch ${prefix}.silhouette.png
    touch ${prefix}.davies_bouldin.png
    touch ${prefix}.calinski_harabasz.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        matplotlib: \$(python3 -c "from importlib.metadata import version; print(version('matplotlib'))")
        numpy: \$(python3 -c "from importlib.metadata import version; print(version('numpy'))")
        pandas: \$(python3 -c "from importlib.metadata import version; print(version('pandas'))")
        scikit-learn: \$(python3 -c "from importlib.metadata import version; print(version('scikit-learn'))")
    END_VERSIONS
    """
}
