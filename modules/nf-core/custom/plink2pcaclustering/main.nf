process CUSTOM_PLINK2PCACLUSTERING {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a3/a37807bdaf3edad30a2b212962b6af381bc10381a80c40efb2bb07f6ee43032f/data' :
        'community.wave.seqera.io/library/numpy_pandas_python_pyyaml_scikit-learn:c500ceb82d3d7606' }"

    input:
    tuple val(meta), path(eigenvec)
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
