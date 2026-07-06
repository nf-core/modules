process RIKER_WGS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/62ce363fc85eaa178522adfc3ddbc4a145d7a4136981e060151886e34e7a55d5/data' :
        'community.wave.seqera.io/library/riker:0.3.0--56fa17ae2be0828f' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)
    path intervals

    output:
    tuple val(meta), path("*.wgs-metrics.txt")  , emit: metrics
    tuple val(meta), path("*.wgs-coverage.txt") , emit: coverage
    tuple val(meta), path("*.wgs-coverage.pdf") , emit: pdf
    tuple val("${task.process}"), val('riker'), eval("riker --version 2>&1 | sed 's/riker //'") , topic: versions, emit: versions_riker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def interval_arg = intervals ? "-L ${intervals}" : ''
    """
    riker wgs \\
        -i ${bam} \\
        -r ${fasta} \\
        -o ${prefix} \\
        ${interval_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.wgs-metrics.txt
    touch ${prefix}.wgs-coverage.txt
    touch ${prefix}.wgs-coverage.pdf
    """
}
