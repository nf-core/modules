// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SPATIALCONTINUUM_CRICALC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a90aa9a3b3a71e63730c7004c9516005e4080233056e7560e28467ea0dbc726b/data'
        : 'community.wave.seqera.io/library/spatial-continuum:1.0.3--ae8c7c8dbbc78ef7'}"

    input:
    tuple val(meta), val(pattern), val(tp53), val(rb1), val(stroma)

    output:
    tuple val(meta), path("*.txt"), emit: report
    tuple val("${task.process}"), val('spatial-continuum'), eval("python3 -c \"import importlib.metadata; print(importlib.metadata.version('spatial-continuum'))\""), topic: versions, emit: versions_spatialcontinuum

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cri-calc \\
        --pattern "${pattern}" \\
        --tp53 "${tp53}" \\
        --rb1 "${rb1}" \\
        --stroma "${stroma}" \\
        ${args} \\
        > ${prefix}_cri_report.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cri_report.txt
    """
}
