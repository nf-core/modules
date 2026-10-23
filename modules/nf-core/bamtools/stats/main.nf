// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c1/c1fde4cd562bddffee839570552212cbc3e01bb253ae8af6954a069106d518ba/data'
        : 'community.wave.seqera.io/library/bamtools:2.5.2--b804fdcee08dfdd6' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('bamtools'), eval("bamtools --version | sed '2!d;s/bamtools //g'"), emit: versions_bamtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamtools \\
        stats \\
        -in $bam \\
        $args \\
        >${prefix}.bam.stats
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam.stats
    """
}
