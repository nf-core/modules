// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BAMALIGNCLEANER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/857ab804ba29aedc4b653a55fe9318a182140140042cc141468a1688bdf65ac9/data'
        : 'community.wave.seqera.io/library/bamaligncleaner:0.2.2--8e64749f97f80978' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('bamaligncleaner'), eval("bamAlignCleaner --version | sed 's/.*version //'"), emit: versions_bamaligncleaner, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bamAlignCleaner \\
        $args \\
        -o ${prefix}.bam \\
        ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    """
}
