// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEBREAK {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/10703319a506f65370b48d639d35334a8b281ccfb9984dee22ae02ae531c3303/data'
        : 'community.wave.seqera.io/library/debreak:1.3--7258ebe5834e7801'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('debreak'), val('1.3'), emit: versions_debreak, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args ?: ''
    def ref_arg = fasta ? "-r ${fasta}" : ""
    """
    debreak \\
        --bam ${bam} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${ref_arg} \\
        ${args}

    mv ${prefix}/debreak.vcf ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
