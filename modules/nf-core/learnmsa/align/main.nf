// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process LEARNMSA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca479ece475f5b115a5e267275ad066fba450c2a4d2527cb6ca7e591a8ff8a7f/data'
        : 'community.wave.seqera.io/library/learnmsa:2.0.16--7965085f02b1399a'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aln"), emit: alignment
    tuple val("${task.process}"), val('learnmsa'), eval("learnMSA -h | sed -nE 's/.*version ([0-9.]+).*/\\1/p'"), emit: versions_learnmsa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    learnMSA \\
        -i ${fasta} \\
        -o "${prefix}.aln" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln
    """
}
