// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/44/4433b2cb821d4bc152f1911ab5750ca0feda7980b5e3d87a3eccbf864d0821ee/data'
        : 'community.wave.seqera.io/library/bioconductor-summarizedexperiment:1.42.0--39b9fc09268fd17d'}"

    input:
    tuple val(meta), path(matrix_files)
    tuple val(meta2), path(rowdata)
    tuple val(meta3), path(coldata)

    output:
    tuple val(meta), path("*.rds")              , emit: rds
    tuple val(meta), path("*.R_sessionInfo.log"), emit: log
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'summarizedexperiment.r'

    stub:
    """
    touch ${meta.id}.SummarizedExperiment.rds
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
