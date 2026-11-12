// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SEQTK_COMP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2dae00435c9dff9202076069bf0d252db26a3e8c708741f17b5fbd053886ff9a/data'
:         'community.wave.seqera.io/library/seqtk:1.4--3fa479d3961a5b8d' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.seqtk_stats.tsv"), emit: seqtk_stats
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | sed -n 's/^Version: //p'"), emit: versions_seqtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqtk comp \\
        ${args} \\
        ${fastx} > ${prefix}.seqtk_stats.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" > ${prefix}.seqtk_stats.tsv
    """
}
