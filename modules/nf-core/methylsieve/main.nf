// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process METHYLSIEVE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b6ddf32365a6fe29b43b89e0ed11ca9dda2d3f73d9a72827c9910a24ef169b2f/data'
        : 'community.wave.seqera.io/library/methylsieve:0.3.0--fabcb86b0e4a4051'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("*.summary.tsv"), emit: summary
    tuple val(meta), path("*.mbias.tsv"), emit: mbias
    tuple val(meta), path("*.conversion-matrix.tsv"), emit: conversion_matrix, optional: true
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val("${task.process}"), val('methylsieve'), eval("methylsieve --version | sed 's/^methylsieve //'"), emit: versions_methylsieve, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    methylsieve \\
        --reference ${fasta} \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --metrics-prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.summary.tsv
    touch ${prefix}.mbias.tsv
    touch ${prefix}.conversion-matrix.tsv
    touch ${prefix}.mbias.pdf
    touch ${prefix}.conversion-matrix.pdf
    """
}
