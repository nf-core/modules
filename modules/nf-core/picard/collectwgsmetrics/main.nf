process PICARD_COLLECTWGSMETRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path intervallist

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val("${task.process}"), val('picard'), eval("picard CollectWgsMetrics --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    def interval = intervallist ? "--INTERVALS ${intervallist}" : ''
    if (!task.memory) {
        log.info('[Picard CollectWgsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectWgsMetrics \\
        ${args} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.CollectWgsMetrics.coverage_metrics \\
        --REFERENCE_SEQUENCE ${fasta} \\
        ${interval}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectWgsMetrics.coverage_metrics
    """
}
