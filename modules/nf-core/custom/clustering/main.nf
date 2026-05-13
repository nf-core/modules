process CUSTOM_CLUSTERING {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'docker://community.wave.seqera.io/library/matplotlib_pandas_python_scikit-learn_pruned:054f91aaa56bd7d5' :
    'community.wave.seqera.io/library/matplotlib_pandas_python_scikit-learn_pruned:054f91aaa56bd7d5' }"

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
    template 'clustering.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clusters.csv
    touch ${prefix}_clustering_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
