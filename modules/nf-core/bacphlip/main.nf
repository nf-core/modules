// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BACPHLIP {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0'
        : 'quay.io/biocontainers/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bacphlip"), emit: bacphlip_results
    tuple val(meta), path("*.hmmsearch.tsv"), emit: hmmsearch_results
    tuple val("${task.process}"), val('bacphlip'), val('0.9.6'), emit: versions_bacphlip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bacphlip \\
        -i ${fasta} \\
        ${args}
    """

    stub:
    """
    touch ${fasta}.bacphlip
    touch ${fasta}.hmmsearch.tsv
    """
}
