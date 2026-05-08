process CUSTOM_CLUSTERMETRICS {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/69/69a6d33f6bd1a901cad8a6914b6ad11a7db6c35005b4ff8604f20f1baba10fc3/data' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_scikit-learn:b7d7028d28dc4084' }"
    input:
    tuple val(meta), path(features), path(clusters)
    val out_prefix

    output:
    tuple val(meta), path("*_metrics.tsv")     , emit: metrics
    tuple val(meta), path("*_k_sweep.csv")     , emit: k_sweep
    tuple val(meta), path("*_selected.json")   , emit: selected
    tuple val(meta), path("*.png")             , emit: plots, optional: true
    path "versions.yml"                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'cluster_metrics.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metrics.tsv
    touch ${prefix}_k_sweep.csv
    touch ${prefix}_selected.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
