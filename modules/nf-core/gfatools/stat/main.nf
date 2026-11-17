// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GFATOOLS_STAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/590e92e087843842a825609903dc3b83bba5077f768aafd48af0f551ddb143e8/data'
        : 'community.wave.seqera.io/library/gfatools:0.5--d580643976f2e334'}"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('gfatools'), eval("gfatools version | sed '1!d; s/.* //'"), topic: versions, emit: versions_gfatools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gfatools \\
        stat \\
        $args \\
        $gfa \\
        > ${prefix}.stats
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    """
}
