// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MOTUS_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/58/580145dd529afb23666d3e79ed2dfffce058ab73b55785a23023c104b74cb52b/data'
        : 'community.wave.seqera.io/library/bwa_motus_metasnv_htslib_pruned:813ab30f1d170975'}"

    input:
    tuple val(meta), path(input)
    path db

    output:
    tuple val(meta), path("*.txt") , emit: txt, optional: true
    tuple val(meta), path("*.biom"), emit: biom, optional: true
    // WARN: Version information not provided by tool on CLI.  Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('motus'), val("3.1.0"), topic: versions, emit: versions_motus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd_input = input.size() > 1 ? "-i ${input.join(',')}" : input.isDirectory() ? "-d ${input}" : "-i ${input}"
    def suffix = task.ext.args?.contains("-B") ? "biom" : "txt"
    """
    motus \\
        merge \\
        -db ${db} \\
        ${cmd_input} \\
        ${args} \\
        -o ${prefix}.${suffix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.args?.contains("-B") ? "biom" : "txt"

    """
    touch ${prefix}.${suffix}
    """

}
