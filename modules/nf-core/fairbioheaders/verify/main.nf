// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FAIRBIOHEADERS_VERIFY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a69700529d3e088166672547fbd632315145414a75ee81a1686e1208038a937/data'
        : 'community.wave.seqera.io/library/fair-bioheaders:0.9.0--71029206d3ce1119'}"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${prefix}.verify.txt"), emit: report
    tuple val("${task.process}"), val('fair-bioheaders'), eval("bioheaders --version"), topic: versions, emit: versions_fairbioheaders

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // A group, not a pipe, so that a failed check fails the task.
    // Nextflow already shell-escapes staged inputs, so only the prefix is quoted.
    """
    {
        printf 'file: %s\\n' ${sequence}
        bioheaders \\
            verify \\
            ${args} \\
            ${sequence}
    } > "${prefix}.verify.txt"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.verify.txt"
    """
}
