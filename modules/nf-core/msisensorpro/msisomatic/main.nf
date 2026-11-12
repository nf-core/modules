// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MSISENSORPRO_MSISOMATIC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/78b9be849e5052505c2c7369b8c080dd0cfddb239c7f48f7a9fef65bc6e88dcd/data'
        : 'community.wave.seqera.io/library/msisensor-pro:1.3.0--31bc4af4e656f312'}"

    input:
    tuple val(meta), path(normal), path(normal_index), path(tumor), path(tumor_index), path(intervals)
    tuple val(meta2), path(fasta)
    path(msisensor_scan)

    output:
    tuple val(meta), path("${prefix}")         , emit: output_report
    tuple val(meta), path("${prefix}_dis")     , emit: output_dis
    tuple val(meta), path("${prefix}_all")     , emit: output_all
    tuple val(meta), path("${prefix}_unstable"), emit: output_unstable
    tuple val("${task.process}"), val('msisensor-pro'), eval("msisensor-pro --version 2>&1 | sed -nE 's/Version:\\s*v//p'") , emit: versions_msisensorpro, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def fasta_cmd     = fasta     ? "-g ${fasta}"       : ""
    def intervals_cmd = intervals ? " -e ${intervals} " : ""

    """
    msisensor-pro \\
        msi \\
        -d ${msisensor_scan} \\
        -n ${normal} \\
        -t ${tumor} \\
        ${fasta_cmd} \\
        -o ${prefix} \\
        -b ${task.cpus} \\
        ${intervals_cmd} \\
        ${args}
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}
    touch ${prefix}_dis
    touch ${prefix}_all
    touch ${prefix}_unstable
    """
}
