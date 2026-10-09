// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SPRING_COMPRESS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f67f27c8cb2d1a149564f1a10f5f2b7a6acfa87ef3d3d27d2d8752dbe95e6acf/data'
        : 'community.wave.seqera.io/library/spring:1.1.1--911a17b4ccfb85ee'}"

    input:
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta), path("*.spring"), emit: spring
    tuple val("${task.process}"), val('spring'), val('1.1.1'), topic: versions, emit: versions_spring

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-i ${fastq1}" : "-i ${fastq1} ${fastq2}"
    """
    spring \\
        -c \\
        -g \\
        -t ${task.cpus} \\
        ${args} \\
        ${input} \\
        -o ${prefix}.spring
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" > ${prefix}.spring
    """
}
