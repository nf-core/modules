process CUSTOM_CLUSTERING {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d4df0b5c9756861b31d9f3db7771da8119c9f9646724bdc1c822672db57af93e/data' :
        'community.wave.seqera.io/library/numpy_pandas_python_scikit-learn:ef102bc9b3784075' }"

    input:
    tuple val(meta), path(eigenvec)
    val algorithm
    val n_clusters
    val dbscan_eps
    val dbscan_min_samples

    output:
    tuple val(meta), path("*_clusters.csv")        , emit: clusters
    tuple val(meta), path("*_clustering_info.json"), emit: info
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'clustering.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clusters.csv
    touch ${prefix}_clustering_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
