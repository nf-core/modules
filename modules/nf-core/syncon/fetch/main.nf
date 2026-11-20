// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SYNCON_FETCH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/be85cb074d7868dc6fbfcdf0e6eabde839e61021665c8c2f856aa96dcbc3e568/data'
        : 'community.wave.seqera.io/library/syncon:1.0--454827d120667dec'}"

    input:
    tuple val(meta), path(sequences), path(decon_txt)

    output:
    tuple val(meta), path("*.fa"), emit: decontaminated_fasta
    // WARN: Incorrect version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('syncon'), val("1.0"), topic: versions, emit: versions_syncon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    syncon fetch \\
        ${args} \\
        -i <(grep KEEP ${decon_txt}) \\
        -o ${prefix}.fa \\
        ${sequences} \\
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.fa
    """
}
