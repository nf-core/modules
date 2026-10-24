// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BAMUTIL_CLIPOVERLAP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91c86394353680e525a569121749574c1340e4b82168af825e4cbed716f9cfc1/data'
        : 'community.wave.seqera.io/library/bamutil:1.0.15--1582863f9eda77c6' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: stats_log, optional: true
    tuple val("${task.process}"), val('bamutil'), eval("bam clipOverlap 2>&1 | grep '^Version:' | sed 's/^Version: //;s/;.*//'"), emit: versions_bamutil, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_clipoverlap"
    """
    bam \\
        clipOverlap \\
        --in ${bam} \\
        --out ${prefix}.bam \\
        --stats \\
        ${args} \\
        2> >( tee ${prefix}.log >&2 ) \\
        | tee ${prefix}.err
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_clipoverlap"
    """
    touch ${prefix}.bam
    touch ${prefix}.log
    """
}
