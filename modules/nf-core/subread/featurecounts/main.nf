// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SUBREAD_FEATURECOUNTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/49d1257914d0c2c6bbf05d04b8bad37f23c20369a269033d9973fca590da89ca/data'
        : 'community.wave.seqera.io/library/subread:2.1.1--c888fa34748dca4c'}"

    input:
    tuple val(meta), path(bams), path(annotation)

    output:
    tuple val(meta), path("*featureCounts.tsv"), emit: counts
    tuple val(meta), path("*featureCounts.tsv.summary"), emit: summary
    tuple val("${task.process}"), val('subread'), eval("featureCounts -v 2>&1 | sed 's/featureCounts v//'"), emit: versions_subread, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-p'

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    }
    else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    """
    featureCounts \\
        ${args} \\
        ${paired_end} \\
        -T ${task.cpus} \\
        -a ${annotation} \\
        -s ${strandedness} \\
        -o ${prefix}.featureCounts.tsv \\
        ${bams.join(' ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.tsv
    touch ${prefix}.featureCounts.tsv.summary
    """
}
