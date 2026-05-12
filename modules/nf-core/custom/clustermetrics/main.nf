process CUSTOM_CLUSTERMETRICS {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'docker://community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:6b81abc92579656a' :
    'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:6b81abc92579656a' }"
    input:
    tuple val(meta), path(features), path(clusters)

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
    touch ${prefix}_elbow.png
    touch ${prefix}_silhouette.png
    touch ${prefix}_davies_bouldin.png
    touch ${prefix}_calinski.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
