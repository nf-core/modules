// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TAUGLYMPHATIC_CALC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db230f0ed6fef58e6960f9a04a15307470749a848c165ea6c2964803d5505471/data'
        : 'community.wave.seqera.io/library/tau-glymphatic-continuum:1.1.2--19d94451d3401cc3'}"

    input:
    val meta

    output:
    tuple val(meta), path("${prefix}.txt") , emit: report, optional: true
    tuple val(meta), path("${prefix}.json"), emit: json  , optional: true
    tuple val("${task.process}"), val('tau-glymphatic-continuum'), eval("python3 -c \"import importlib.metadata; print(importlib.metadata.version('tau-glymphatic-continuum'))\""), topic: versions, emit: versions_tauglymphatic

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ext = args.contains('--json') ? 'json' : 'txt'
    """
    tau-calc \\
        ${args} \\
        > ${prefix}.${ext}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ext = args.contains('--json') ? 'json' : 'txt'
    """
    touch ${prefix}.${ext}
    """
}
