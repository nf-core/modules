process SCDBLFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'PLACEHOLDER_SINGULARITY'
        : 'PLACEHOLDER_DOCKER'}"

    input:
    tuple val(meta), path(h5ad), val(dbr), val(batch_col)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.csv") , emit: predictions
    path "versions.yml"            , emit: versions_scdblfinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_scdblfinder"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template 'scdblfinder.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_scdblfinder"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.h5ad
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        scdblfinder: \$(Rscript -e "suppressMessages(library(scDblFinder)); cat(as.character(packageVersion('scDblFinder')))")
    END_VERSIONS
    """
}
