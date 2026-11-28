// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process INFLAMMATORYONCOPHYSICS_CALC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/69/69015543f83041aebf808b0a6c411cfedba10764ec0ff871a2d6c3d5e8a93922/data'
        : 'community.wave.seqera.io/library/inflammatory-onco-physics:1.1.2--b2972c171ba3d1ee'}"

    input:
    tuple val(meta), val(stage), val(grade), val(nlr), val(crp), val(alb)

    output:
    tuple val(meta), path("*.txt"), emit: report
    tuple val("${task.process}"), val('inflammatory-onco-physics'), eval("python3 -c \"import importlib.metadata; print(importlib.metadata.version('inflammatory-onco-physics'))\""), topic: versions, emit: versions_inflammatoryoncophysics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    oncophysics-calc \\
        --stage ${stage} \\
        --grade ${grade} \\
        --nlr ${nlr} \\
        --crp ${crp} \\
        --alb ${alb} \\
        ${args} \\
        > ${prefix}_htbi_report.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_htbi_report.txt
    """
}
