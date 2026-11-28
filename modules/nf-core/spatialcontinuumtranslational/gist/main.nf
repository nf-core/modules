// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SPATIALCONTINUUMTRANSLATIONAL_GIST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a50a8bbb2134cf3a6b6aa7214ba42989fda33ed2804a1bf0bcf2d24dd5f62cce/data'
        : 'community.wave.seqera.io/library/spatial-continuum-translational:0.1.53--ee1f749133648d82'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${prefix}_gist.json"), emit: json, optional: true
    tuple val(meta), path("${prefix}_gist.tsv") , emit: tsv,  optional: true
    tuple val("${task.process}"), val('spatial-continuum-translational'), eval("python3 -c \"import importlib.metadata; print(importlib.metadata.version('spatial-continuum-translational'))\""), topic: versions, emit: versions_spatialcontinuumtranslational

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def format_match = (args =~ /(?:--format|-f)[=\s]+['"]?([a-zA-Z0-9]+)['"]?/)
    def format = format_match ? format_match[0][1].toLowerCase() : 'json'
    """
    spatial-continuum-gist \\
        --input "${input}" \\
        --output "${prefix}_gist.${format}" \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def format_match = (args =~ /(?:--format|-f)[=\s]+['"]?([a-zA-Z0-9]+)['"]?/)
    def format = format_match ? format_match[0][1].toLowerCase() : 'json'
    """
    touch "${prefix}_gist.${format}"
    """
}
