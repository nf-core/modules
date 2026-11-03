// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ABRICATE_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b75992c3e2d82899245f0e5839007cc1dc053cfa72d6018a66729e2ba540af75/data'
        : 'community.wave.seqera.io/library/abricate:1.0.1--cb0de91fb575483f'}"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("*.txt"), emit: report
    tuple val("${task.process}"), val('abricate'), eval("abricate --version 2>&1 | sed 's/^.*abricate //'"), emit: versions_abricate_summary, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    abricate \\
        --summary \\
        ${reports} > ${prefix}.txt
    """
}
