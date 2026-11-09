// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HAMRONIZATION_RESFINDER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3ffa68bd667eec67b6efbf2bf8aeeafd7f69228d537959be178249d021018817/data':
        'community.wave.seqera.io/library/hamronization:1.3.1--86ff1149c2d6acff' }"

    input:
    tuple val(meta), path(report)
    val format

    output:
    tuple val(meta), path("*.json"), optional: true, emit: json
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    tuple val("${task.process}"), val('hamronization'), eval("hamronize --version 2>&1 | sed 's/hamronize //'"), topic: versions, emit: versions_hamronization

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hamronize \\
        resfinder \\
        ${report} \\
        ${args} \\
        --format ${format} \\
        > ${prefix}.${format}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${format}
    """
}
