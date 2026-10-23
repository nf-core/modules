// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BEDOPS_CONVERT2BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46dc39b50ab3b00b757de451db624d1fdec87ee14a0cd5234598b71f200d92e0/data'
        : 'community.wave.seqera.io/library/bedops:2.4.42--1f31d470625cb9a3' }"

    input:
    tuple val(meta), path(in_file)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val("bedops"), eval('convert2bed --version | sed -n "s/.*version: *\\([^ ]*\\).*/\\1/p"'), emit: versions_bedops, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = in_file.getExtension()
    """
    convert2bed \\
        $args \\
        -i $format \\
        < $in_file \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
