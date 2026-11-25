// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEACON_INDEXINTERSECT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f1270fe75323ac6d213bf66d72b9f292c5ff31f7f5ee12300d2bed2b858f5ab/data'
        : 'community.wave.seqera.io/library/deacon:0.20.0--09d7e8e5c3f9dc88'}"

    input:
    tuple val(meta), path(indices)  // two or more deacon .idx index files to intersect

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.intersect"

    """
    deacon \\
        index \\
        intersect \\
        ${args} \\
        ${indices.join(' ')} \\
        > ${prefix}.idx
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.intersect"
    """
    touch ${prefix}.idx
    """
}
