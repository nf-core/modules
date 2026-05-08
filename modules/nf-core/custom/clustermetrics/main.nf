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
    def prefix = task.ext.prefix ?: out_prefix ?: "${meta.id}"

    """
    python3 ${moduleDir}/templates/cluster_metrics.py \\
        --features ${features} \\
        --clusters ${clusters} \\
        --out-k-sweep ${prefix}_k_sweep.csv \\
        --out-selected ${prefix}_selected.json \\
        --out-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
