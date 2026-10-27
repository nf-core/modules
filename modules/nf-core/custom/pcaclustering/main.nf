// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_PCACLUSTERING {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7954a01492f2eda956d7646dfdc88f3146e2a65cfd2c4cb22f935938408811cd/data'
        : 'community.wave.seqera.io/library/numpy_pandas_python_pyyaml_scikit-learn:b745a964522e2e96'}"

    input:
    tuple val(meta), path(features)
    val algorithm
    val n_clusters
    val dbscan_eps
    val dbscan_min_samples

    output:
    tuple val(meta), path("*.clusters.csv")        , emit: clusters
    tuple val(meta), path("*.clustering_info.json"), emit: info
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'clustering.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clusters.csv
    touch ${prefix}.clustering_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        numpy: \$(python3 -c "from importlib.metadata import version; print(version('numpy'))")
        pandas: \$(python3 -c "from importlib.metadata import version; print(version('pandas'))")
        scikit-learn: \$(python3 -c "from importlib.metadata import version; print(version('scikit-learn'))")
    END_VERSIONS
    """
}
