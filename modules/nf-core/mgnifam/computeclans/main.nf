// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MGNIFAM_COMPUTECLANS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07d734ffff975d82e7381943964ea0b875606d69d7f3b3c35d1488833fff960b/data'
        : 'community.wave.seqera.io/library/mgnifam:4.1.0--818ca9d3b99bfbca'}"

    input:
    tuple val(meta), path(members), path(fasta_file)

    output:
    tuple val(meta), path("${prefix}/${prefix}_clans.csv")         , emit: clans
    tuple val(meta), path("${prefix}/${prefix}_clan_edges.csv")    , emit: clan_edges
    tuple val(meta), path("${prefix}/${prefix}_clans.log")         , emit: log
    tuple val(meta), path("${prefix}/${prefix}_clans_summary.json"), emit: summary
    tuple val("${task.process}"), val('mgnifam'), eval("mgnifam --version 2>&1"), topic: versions, emit: versions_mgnifam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def fasta = fasta_file ? "--fasta_file \"${fasta_file}\"" : ''
    """
    mgnifam compute_clans \\
        --members ${members} \\
        --output_dir "${prefix}" \\
        --cpus ${task.cpus} \\
        --chunk_id "${prefix}" \\
        ${fasta} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}
    touch ${prefix}/${prefix}_clans.csv
    touch ${prefix}/${prefix}_clan_edges.csv
    touch ${prefix}/${prefix}_clans.log
    touch ${prefix}/${prefix}_clans_summary.json
    """
}
