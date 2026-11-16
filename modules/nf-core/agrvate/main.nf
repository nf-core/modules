// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AGRVATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96dcafb764303f639281ee284658f92cb2abced5fb6285bdfcc965763ae98fef/data'
        : 'community.wave.seqera.io/library/agrvate:1.0.2--a11be0fee842dcd2'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}-results/${fasta.baseName}-summary.tab"), emit: summary
    path "${fasta.baseName}-results"                                                , emit: results_dir
    tuple val("${task.process}"), val('agrvate'), eval("agrvate --version | sed 's/[^0-9.]//g'"), emit: versions_agrvate, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    agrvate \\
        ${args} \\
        -i ${fasta}
    """

    stub:

    """
    mkdir ${fasta.baseName}-results
    touch ${fasta.baseName}-results/${fasta.baseName}-summary.tab
    """
}
