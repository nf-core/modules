process SCDBLFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b7da7136d1cb87691c20da1e7b64b2b9bbd102baee2fbc803aedb7c5f634e7f/data'
        : 'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-biocparallel_bioconductor-rhdf5_bioconductor-scdblfinder_bioconductor-singlecellexperiment:b97c9b42c2a89cba'}"

    input:
    tuple val(meta), path(h5ad), val(dbr), val(batch_col)

    output:
    tuple val(meta), path("*.h5ad")            , emit: h5ad
    tuple val(meta), path("*.predictions.csv") , emit: predictions
    tuple val(meta), path("*.scores.csv")      , emit: scores
    path "versions.yml"                        , emit: versions_scdblfinder, topic: versions

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
    touch ${prefix}.predictions.csv
    touch ${prefix}.scores.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        scdblfinder: \$(Rscript -e "suppressMessages(library(scDblFinder)); cat(as.character(packageVersion('scDblFinder')))")
    END_VERSIONS
    """
}
