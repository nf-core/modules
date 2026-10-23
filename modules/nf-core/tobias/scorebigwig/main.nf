// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TOBIAS_SCOREBIGWIG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb99dc2a6b4da1bed5ad4d66f889e3830d7682d9e727df26337434d784868353/data'
        : 'community.wave.seqera.io/library/tobias_python:33fc15661f1e50e6'}"

    input:
    tuple val(meta), path(signal), path(regions)

    output:
    tuple val(meta), path("*_footprints.bw"), emit: footprints
    tuple val("${task.process}"), val('tobias'), eval('TOBIAS --version'), topic: versions, emit: versions_tobias

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p .matplotlib
    export MPLCONFIGDIR="\${PWD}/.matplotlib"

    TOBIAS ScoreBigwig \\
        --signal $signal \\
        --regions $regions \\
        --output ${prefix}_footprints.bw \\
        --cores ${task.cpus} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_footprints.bw
    """
}
