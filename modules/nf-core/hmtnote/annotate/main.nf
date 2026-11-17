// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HMTNOTE_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a900304b834f8965ad77c0976311b69f97be9f96f03cb04e83082ca0d265bdd7/data'
        : 'community.wave.seqera.io/library/hmtnote_python:6774601e4144aa5d'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_annotated.vcf"), emit: vcf
    tuple val("${task.process}"), val('hmtnote'), eval("hmtnote --version 2>&1 | sed 's/.*version //'"), emit: versions_hmtnote, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated.

Reason:
This tool is no longer maintained by its author and, as its database hosting service has
been discontinued, it can no longer work with conda.
"""
    assert false: deprecation_message
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    hmtnote \\
        annotate \\
        ${vcf} \\
        ${prefix}_annotated.vcf \\
        ${args}
    """

    stub:
        def deprecation_message = """
WARNING: This module has been deprecated.

Reason:
This tool is no longer maintained by its author and, as its database hosting service has
been discontinued, it can no longer work with conda.
"""
    assert false: deprecation_message
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.vcf
    """
}
