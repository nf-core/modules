// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process METAPHLAN_MERGEMETAPHLANTABLES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/59f3aa387903d5b87b6e2d3c0172f462f701c662e73cb2e26066d4ec6fd3e8c6/data'
        : 'community.wave.seqera.io/library/metaphlan_numpy_python:59796ca881a8bb02'}"

    input:
    tuple val(meta), path(profiles)

    output:
    tuple val(meta), path("${prefix}.txt") , emit: txt
    tuple val("${task.process}"), val('metaphlan'), eval("metaphlan --version 2>&1 | cut -d ' ' -f 3"), emit: versions_metaphlan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_metaphlan_tables.py \\
        $args \\
        -o ${prefix}.txt \\
        ${profiles}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    """
}
