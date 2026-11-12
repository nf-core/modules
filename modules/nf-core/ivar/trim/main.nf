// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf3ee0025f72df3ca5cabfae265b3623dcd9f9239f993ea258bf5e7846e9a4e4/data'
        : 'community.wave.seqera.io/library/ivar:1.4.4--60def7bdef0d4d2f'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path('*.log'), emit: log
    tuple val("${task.process}"), val('ivar'), eval("ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p'"), emit: versions_ivar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar trim \\
        ${args} \\
        -i ${bam} \\
        -b ${bed} \\
        -p ${prefix} \\
        > ${prefix}.ivar.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ivar.log
    touch ${prefix}.bam
    """
}
