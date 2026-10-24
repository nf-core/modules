// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FGUMI_CLIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/74/74c00516b0152fbf7f721f7b8983024c652639d47c5e225c0ad0f6ca00773788/data'
:         'community.wave.seqera.io/library/fgumi:0.8.0--4fd382d441fde4a1' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics, optional: true
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_clipped"
    // fgumi clip cannot combine --metrics with --threads, so only request threads when metrics are not
    def threads = args.contains('--metrics') ? '' : "--threads ${task.cpus}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    fgumi \\
        clip \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --reference ${fasta} \\
        ${threads} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_clipped"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    """
}
