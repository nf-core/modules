// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ATAQV_MKARV {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7063be8c7419222c75482548d138078837657f628d149e1a0b3e6a9080799de/data'
        : 'community.wave.seqera.io/library/ataqv:1.3.1--13b9c1b1eb8ada3f'}"

    input:
    path "jsons/*"

    output:
    path "html"        , emit: html
    tuple val("${task.process}"), val('ataqv'), eval("ataqv --version 2>&1 || true"), emit: versions_ataqv, topic: versions
    // tuple val("${task.process}"), val('mkarv'), eval('mkarv --version'), emit: versions_mkarv, topic: versions //Use this when version string has been fixed
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkarv \\
        ${args} \\
        --concurrency ${task.cpus} \\
        --force \\
        ./html/ \\
        jsons/*
    """

    stub:
    """
    mkdir -p html
    touch html/index.html
    """
}
