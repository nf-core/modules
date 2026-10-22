process CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/celltypist:1.7.1--pyhdfd78af_2'
        : 'quay.io/biocontainers/celltypist:1.7.1--pyhdfd78af_2'}"

    input:
    tuple val(meta), path(h5ad)
    tuple val(meta2), path(model)

    output:
    tuple val(meta), path("*predicted_labels.csv")  , emit: predicted_labels
    tuple val(meta), path("*decision_matrix.csv")   , emit: decision_matrix
    tuple val(meta), path("*probability_matrix.csv"), emit: probability_matrix
    tuple val("${task.process}"), val('celltypist'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"celltypist\"))'"), topic: versions, emit: versions_celltypist

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Home is read-only in some containers, so keep CellTypist's cache and plotting files in the work directory
    export CELLTYPIST_FOLDER=./celltypist_cache
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    celltypist \\
        ${args} \\
        --indata ${h5ad} \\
        --model ${model} \\
        --outdir . \\
        --prefix ${prefix}_
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_labels.csv
    touch ${prefix}_decision_matrix.csv
    touch ${prefix}_probability_matrix.csv
    """
}
