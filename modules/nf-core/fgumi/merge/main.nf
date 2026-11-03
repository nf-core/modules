// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FGUMI_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8b/8b1b1f653fc1fef5eb7d31f72d7cc60e29ff0eec5552802523ce97b1b5d0f920/data'
        : 'community.wave.seqera.io/library/fgumi:0.9.0--1f0a2f2a31a251a3'}"

    input:
    tuple val(meta), path(bams, stageAs: "?/*")

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgumi \\
        merge \\
        --output ${prefix}.bam \\
        --threads ${task.cpus} \\
        ${args} \\
        ${bams}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
