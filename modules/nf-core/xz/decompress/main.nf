// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process XZ_DECOMPRESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db54d725f04a667223c103459fb70612a72f9766bfd74b26f2d566e1b190a65d/data'
:         'community.wave.seqera.io/library/xz:5.2.6--5bcc1ddb5aafa435' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$decompressed_file"), emit: file
    tuple val("${task.process}"), val('xz'), eval('xz --version | sed \'1!d;s/\\([^0-9.]*\\)//g\''), topic: versions, emit: versions_xz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    decompressed_file = archive.toString().replaceAll(".xz\$", "")
    """
    # Note 1: needs --stdout for xz --decompress to avoid two issues:
    #   1. xz: ${archive}: Is a symbolic link, skipping
    #   2. xz: ${archive}: Cannot set the file group: Operation not permitted
    # Note 2: using several threads in xz --decompress will only work on files that contain multiple blocks with size
    # information in block headers.  All files compressed in multi-threaded mode meet this condition.
    xz -T ${task.cpus} --decompress --stdout ${args} ${archive} > ${decompressed_file}
    """

    stub:
    def args = task.ext.args ?: ''
    decompressed_file = archive.toString().replaceAll(".xz\$", "")
    """
    echo "${args}"
    touch "${decompressed_file}"
    """
}
