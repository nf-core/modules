// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ASSEMBLYSCAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f8b13fe9789963debbbbe0e934aea5e96cfef77277299e3152726738fb81444a/data'
        : 'community.wave.seqera.io/library/assembly-scan:1.0.0--fdab3a3c6ba9bbd0'}"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.{json,tsv}"),     emit: report
    tuple val("${task.process}"), val('assemblyscan'), eval("assembly-scan --version 2>&1 | sed 's/^.*assembly-scan //; s/Using.*\$//'"), topic: versions, emit: versions_assemblyscan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def file_format = (task.ext.args?.contains("--json")) ? "json" : "tsv"
    """
    assembly-scan \\
        ${args} \\
        ${assembly} > ${prefix}.${file_format}
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def file_format = (task.ext.args?.contains("--json")) ? "json" : "tsv"
    """
    touch ${prefix}.${file_format}
    """
}
