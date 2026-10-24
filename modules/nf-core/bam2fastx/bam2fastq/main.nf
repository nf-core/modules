// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BAM2FASTX_BAM2FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/8792febff8eb10b2c43ffdea0bf5a45222b82fb8b7088eb411303d2b5ef5b8fc/data'
        : 'community.wave.seqera.io/library/bam2fastx:1.3.1--1181601d78cb6c32' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('bam2fastx'), eval("bam2fastq --version 2>&1) | sed 's/^.*bam2fastq //'"), emit: versions_bam2fastx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/pbtk/bam2fastq

Reason:
This module is no longer fit for purpose because bam2fastx has been deprecated by PacificBiosciences

"""
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert false: deprecation_message
    """
    bam2fastq \\
        $args \\
        -o ${prefix} \\
        $bam
    """
}
