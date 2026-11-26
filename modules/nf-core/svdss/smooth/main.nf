// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SVDSS_SMOOTH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/da/daa86b771119a67c6dfc0892e58015effe34df2c82754aad558d53253b074600/data'
        : 'community.wave.seqera.io/library/samtools_svdss:a3b4ff04ffb34c8c'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.bam")    , emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai
    tuple val("${task.process}"), val('svdss'), eval("SVDSS --version 2>&1 | sed 's/SVDSS, v//'"), emit: versions_svdss, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.smoothed"
    if ("${bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    SVDSS smooth \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        --bam ${bam} \\
        ${args} \\
        > ${prefix}.bam

    samtools index \\
        -@ ${task.cpus} \\
        ${prefix}.bam
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.smoothed"
    if ("${bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo ${args}

    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    """
}
