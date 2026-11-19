// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SVDSS_SEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8357bb6a78ac936810ecc7528fa4b53a8deea7e7d8abb31afb8ed7e14a05fdc3/data'
        : 'community.wave.seqera.io/library/svdss:2.1.1--5e16609dd2159d4e'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${prefix}.sfs.txt"), emit: sfs
    tuple val("${task.process}"), val('svdss'), eval("SVDSS --version 2>&1 | sed 's/SVDSS, v//'"), emit: versions_svdss, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_flag = reads.name.endsWith('.bam') ? '--bam' : '--fastx'
    """
    SVDSS search \\
        --threads ${task.cpus} \\
        --index ${index} \\
        ${input_flag} ${reads} \\
        ${args} \\
        > ${prefix}.sfs.txt
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.sfs.txt
    """
}
