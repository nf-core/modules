// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SYNCON_BUILD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/be85cb074d7868dc6fbfcdf0e6eabde839e61021665c8c2f856aa96dcbc3e568/data'
        : 'community.wave.seqera.io/library/syncon:1.0--454827d120667dec'}"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.1khash"), emit: onekhash
    tuple val(meta), path("*.1syncseq"), emit: onesyncseq
    // WARN: Incorrect version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('syncon'), val("1.0"), topic: versions, emit: versions_syncon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    syncon build \\
        ${args} \\
        -T ${task.cpus} \\
        -o ${prefix} \\
        ${sequences}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.1khash
    touch ${prefix}.1syncseq
    """
}
