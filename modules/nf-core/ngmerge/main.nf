// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process NGMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c1/c13e2a4a1cfe3fb14a0db0c04f099959036f1cf0de30e1a53193e4ef2c7f2897/data'
        : 'community.wave.seqera.io/library/ngmerge:0.3--b244ad1892587927'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.fq.gz"), emit: merged_reads
    tuple val(meta), path("*_1.fastq.gz")  , emit: unstitched_read1
    tuple val(meta), path("*_2.fastq.gz")  , emit: unstitched_read2
    tuple val("${task.process}"), val('ngmerge'), eval("NGmerge --version 2>&1 | sed -n '1s/^NGmerge, version //p'"), emit: versions_ngmerge, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    NGmerge \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${prefix}.merged.fq.gz \\
        -f ${prefix}_unstitched \\
        -z \\
        -n $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.merged.fq.gz
    """
}
