// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TRUVARI_SEGMENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bcb554578972eea4ca66731f80a41017d591edf30d34f7b3c4a6ef78b0ef4ff3/data'
        : 'community.wave.seqera.io/library/truvari:5.4.0--a783dff809870f2f'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('truvari'), eval("truvari version | sed 's/Truvari v//'"), topic: versions, emit: versions_truvari

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    truvari \\
        segment \\
        -o ${prefix}.vcf \\
        $args \\
        $vcf

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    """
}
