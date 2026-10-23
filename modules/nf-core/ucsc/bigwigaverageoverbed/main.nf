// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process UCSC_BIGWIGAVERAGEOVERBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96d5dedcf4cb3e65e3badb2bd4dabb72beb7fafe919febfbe29a7807a6e60413/data'
:         'community.wave.seqera.io/library/ucsc-bigwigaverageoverbed:482--09b4685dee7c8b33' }"

    input:
    tuple val(meta), path(bed)
    path bigwig

    output:
    tuple val(meta), path("*.tab"), emit: tab
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // BUG: bigWigAverageOverBed cannot handle ensembl seqlevels style
    """
    bigWigAverageOverBed \\
        $args \\
        $bigwig \\
        $bed \\
        ${prefix}.tab
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tab

    """
}
