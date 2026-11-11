// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process COGAPS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b766b163b9a8410332bc4691ea5d738ad2a718dda806e10cc332d5e10fee0a96/data'
        : 'community.wave.seqera.io/library/bioconductor-cogaps:3.33.0--9d268f422491d706'}"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.cogaps.rds")            , emit: rds
    tuple val(meta), path("*.feature_loadings.csv")  , emit: feature_loadings
    tuple val(meta), path("*.sample_factors.csv")    , emit: sample_factors
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'cogaps.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cogaps.rds
    touch ${prefix}.feature_loadings.csv
    touch ${prefix}.sample_factors.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cogaps: \$(Rscript -e 'cat(as.character(packageVersion("CoGAPS")))')
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """
}
