// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_CLUSTERMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b3cc55d016630d9aab5aebf922c59fdc8ef18783ad4572c79ff6b7e4ee2b5964/data'
        : 'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:98f5772f09be086e'}"

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
    touch ${prefix}.k_distance.png

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
