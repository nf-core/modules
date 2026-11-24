// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DUPBLASTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5a7fd2e4aa8b960accf56ecd8c7b66c8e53430188f95a621c95eb19a8c04915e/data'
        : 'community.wave.seqera.io/library/dupblaster:0.3.1--31564d2e7161cbce'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("*.duplicate-metrics.tsv"), emit: metrics
    tuple val(meta), path("*.sequencing-units.tsv"), emit: sequencing_units, optional: true
    tuple val(meta), path("*.duplication-sampled.tsv"), emit: sampled
    tuple val(meta), path("*.duplication-spectrum.tsv"), emit: spectrum, optional: true
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val("${task.process}"), val('dupblaster'), eval("dupblaster --version | sed 's/^dupblaster //'"), emit: versions_dupblaster, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    dupblaster \\
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
    touch ${prefix}.duplicate-metrics.tsv
    touch ${prefix}.sequencing-units.tsv
    touch ${prefix}.duplication-sampled.tsv
    touch ${prefix}.duplication-sampled.pdf
    """
}
