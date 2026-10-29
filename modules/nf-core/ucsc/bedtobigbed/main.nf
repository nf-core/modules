// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39442202396dc06815b4e27d9e2b6af7cb0f46382540da46ebe2304ee15a1f38/data'
:         'community.wave.seqera.io/library/ucsc-bedtobigbed:482--f9a75c5237cbfded' }"

    input:
    tuple val(meta), path(bed)
    path  sizes
    path  autosql

    output:
    tuple val(meta), path("*.bigBed"), emit: bigbed
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def as_option = autosql ? "-as=${autosql}" : ""
    """
    bedToBigBed \\
        $bed \\
        $sizes \\
        $as_option \\
        $args \\
        ${prefix}.bigBed

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bigBed
    """
}
