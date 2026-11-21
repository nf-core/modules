// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FASTCOV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65e3a2366e265fb3fe24cd3082d005be2ab22cde7247a2be651d84a6eaeb07f3/data'
        : 'community.wave.seqera.io/library/fastcov:0.1.3--e3e78d7eeb73e45a'}"

    input:
    tuple val(meta), path(bam), path(index)
    val(file_ext)

    output:
    tuple val(meta), path("${prefix}.${file_ext}")                                                                           , emit: coverage_plot
    tuple val(meta), path("log.txt")                                                                                         , emit: log
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('fastcov'), val('0.1.3')                                                               , emit: versions_fastcov, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    file_ext = file_ext.isEmpty() ? 'png' : file_ext
    """
    fastcov.py \\
        -o ${prefix}.${file_ext} \\
        ${args} \\
        ${bam} > log.txt 2>&1
    """

    stub:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    file_ext = file_ext.isEmpty() ? 'png' : file_ext

    """
    echo $args

    touch ${prefix}.${file_ext}
    touch log.txt
    """
}
