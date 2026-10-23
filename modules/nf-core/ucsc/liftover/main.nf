// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process UCSC_LIFTOVER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15744931def04387cf095376c9d8c322bf21ca542d64946a7258c15b40b12b82/data'
:         'community.wave.seqera.io/library/ucsc-liftover:482--f8639164c0f3b5bf' }"

    input:
    tuple val(meta), path(bed)
    path(chain)

    output:
    tuple val(meta), path("*.lifted.bed")  , emit: lifted
    tuple val(meta), path("*.unlifted.bed"), emit: unlifted
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    liftOver \\
        $args \
        $bed \\
        $chain \\
        ${prefix}.lifted.bed \\
        ${prefix}.unlifted.bed
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lifted.bed
    touch ${prefix}.unlifted.bed
    """

}
