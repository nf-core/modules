process FGUMI_GROUP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a62b457c53300603da026225f95b4db04d1c9f8ba7f734787818fc105d51323/data'
        : 'community.wave.seqera.io/library/fgumi:0.4.0--1fb5dc6de05ce63b'}"

    input:
    tuple val(meta), path(bam)
    val strategy

    output:
    tuple val(meta), path("*.bam")                      , emit: bam
    tuple val(meta), path("*.family_size_histogram.txt"), emit: histogram
    tuple val(meta), path("*.grouping_metrics.txt")     , emit: metrics
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_umi-grouped"

    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    fgumi group \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --strategy ${strategy} \\
        --family-size-histogram ${prefix}.family_size_histogram.txt \\
        --grouping-metrics ${prefix}.grouping_metrics.txt \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_umi-grouped"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.family_size_histogram.txt
    touch ${prefix}.grouping_metrics.txt
    """
}
