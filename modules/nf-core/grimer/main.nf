// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GRIMER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3fa7bb19679b70a3d63ba7724d9357d557249b8ef3d2e2bcc34cf1b654f028c9/data'
        : 'community.wave.seqera.io/library/grimer:1.1.0--3eb74a39b6c2707f'}"

    input:
    tuple val(meta), path(input_table)
    path sample_metadata
    path config

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val("${task.process}"), val('grimer'), eval('grimer --version 2>&1 | tail -1 | sed "s/grimer //"'), emit: versions_grimer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def m_arg  = sample_metadata ? "-m ${sample_metadata}" : ''
    def c_arg  = config          ? "-c ${config}"   : ''
    """
    grimer \\
        -i ${input_table} \\
        ${m_arg} \\
        ${c_arg} \\
        -o ${prefix}.html \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    """
}
