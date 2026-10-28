process CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db26a106d81b1b9dee9d3ce7223c13019a3b8e0ee6797b6197c402637f704017/data'
        : 'community.wave.seqera.io/library/celltypist:1.7.1--81458d01ea674d65'}"

    input:
    tuple val(meta), path(h5ad), path(over_clustering), path(gene_file), path(cell_file)
    tuple val(meta2), path(model)

    output:
    tuple val(meta), path("${prefix}_predicted_labels.csv")  , optional: true, emit: predicted_labels
    tuple val(meta), path("${prefix}_decision_matrix.csv")   , optional: true, emit: decision_matrix
    tuple val(meta), path("${prefix}_probability_matrix.csv"), optional: true, emit: probability_matrix
    tuple val(meta), path("${prefix}_annotation_result.xlsx"), optional: true, emit: xlsx
    tuple val(meta), path("*.pdf")                           , optional: true, emit: plots
    tuple val("${task.process}"), val('celltypist'), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"celltypist\"))'"), topic: versions, emit: versions_celltypist

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def over_clustering_arg = over_clustering ? "--over-clustering ${over_clustering}" : ''
    def gene_file_arg = gene_file ? "--gene-file ${gene_file}" : ''
    def cell_file_arg = cell_file ? "--cell-file ${cell_file}" : ''
    """
    # Home is read-only in some containers, so keep CellTypist's cache and plotting files in the work directory
    export CELLTYPIST_FOLDER=./celltypist_cache
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    celltypist \\
        ${args} \\
        ${over_clustering_arg} \\
        ${gene_file_arg} \\
        ${cell_file_arg} \\
        --indata ${h5ad} \\
        --model ${model} \\
        --outdir . \\
        --prefix ${prefix}_
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tables = args.contains('--xlsx') ? "${prefix}_annotation_result.xlsx" : "${prefix}_predicted_labels.csv ${prefix}_decision_matrix.csv ${prefix}_probability_matrix.csv"
    def plots = args.contains('--plot-results') ? "${prefix}_predicted_labels.pdf" : ''
    """
    touch ${tables} ${plots}
    """
}
