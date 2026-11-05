// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FGUMI_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a01cc9cc4b3bf65b546d5feea63f099a68542b3f4516d4e98e923b13cdba74b/data'
        : 'community.wave.seqera.io/library/fgumi_gzip:f6d152fa9a4b3fcd'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi
    tuple val("${task.process}"), val('gzip'), eval("gzip --version 2>&1 | sed '1!d;s/gzip //'"), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgumi \\
        fastq \\
        --input ${bam} \\
        --threads ${task.cpus} \\
        ${args} \\
        | gzip -c ${args2} > ${prefix}.fastq.gz
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    echo "" | gzip > ${prefix}.fastq.gz
    """
}
