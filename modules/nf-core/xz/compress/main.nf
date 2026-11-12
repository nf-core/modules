// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process XZ_COMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db54d725f04a667223c103459fb70612a72f9766bfd74b26f2d566e1b190a65d/data'
:         'community.wave.seqera.io/library/xz:5.2.6--5bcc1ddb5aafa435' }"

    input:
    tuple val(meta), path(raw_file)

    output:
    tuple val(meta), path("$archive"), emit: archive
    tuple val("${task.process}"), val('xz'), eval('xz --version | sed \'1!d;s/\\([^0-9.]*\\)//g\''), topic: versions, emit: versions_xz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    archive = raw_file.toString() + ".xz"
    """
    # needs --stdout for xz to avoid the following issue:
    # xz: ${raw_file}: Is a symbolic link, skipping
    xz -T $task.cpus --stdout ${args} ${raw_file} > ${archive}
    """

    stub:
    def args = task.ext.args ?: ''
    archive = raw_file.toString() + ".xz"
    """
    echo "${args}"
    touch "${archive}"
    """
}
