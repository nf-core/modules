process FGUMI_DOWNSAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a62b457c53300603da026225f95b4db04d1c9f8ba7f734787818fc105d51323/data':
        'community.wave.seqera.io/library/fgumi:0.4.0--1fb5dc6de05ce63b' }"

    input:
    tuple val(meta), path(bam)
    val fraction
    val keep_rejected

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.rejects.bam"), emit: rejects, optional: true
    tuple val(meta), path("${prefix}.histogram_kept.txt"), emit: histogram_kept
    tuple val(meta), path("${prefix}.histogram_rejected.txt"), emit: histogram_rejected
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_downsampled"
    def rejects_command = keep_rejected ? "--rejects ${prefix}.rejects.bam" : ''
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    fgumi \\
        downsample \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --fraction ${fraction} \\
        --histogram-kept ${prefix}.histogram_kept.txt \\
        --histogram-rejected ${prefix}.histogram_rejected.txt \\
        ${rejects_command} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_downsampled"
    def rejects_command = keep_rejected ? "touch ${prefix}.rejects.bam" : ''
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.histogram_kept.txt
    touch ${prefix}.histogram_rejected.txt
    ${rejects_command}
    """
}
