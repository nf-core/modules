// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GLIMPSE_CHUNK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f722eb1dd56d995b28e1ec47f7808450bf6477e0b4ec28ec5f294c1a66fd479e/data'
        : 'community.wave.seqera.io/library/glimpse-bio:1.1.1--26d63752054f98c5'}"

    input:
    tuple val(meta), path(input), path(input_index), val(region)

    output:
    tuple val(meta), path("*.txt"), emit: chunk_chr
    tuple val("${task.process}"), val('glimpse'), eval("GLIMPSE_chunk --help | sed -n '/Version/s/.*: //p'"), topic: versions, emit: versions_glimpse

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""

    """
    GLIMPSE_chunk \\
        ${args} \\
        --input ${input} \\
        --region ${region} \\
        --thread ${task.cpus} \\
        --output ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt
    """
}
