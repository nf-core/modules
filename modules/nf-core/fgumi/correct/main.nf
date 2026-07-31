process FGUMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64e8f594b6f0dd879bc5abbe4ca70b6b761e1920e407d9e1c7d27b89004aac34/data':
        'community.wave.seqera.io/library/fgumi:0.5.0--a2d14bf52f73eaef' }"

    input:
    tuple val(meta), path(bam)
    path umis
    val min_distance
    val keep_rejected

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.rejects.bam"), emit: rejects, optional: true
    tuple val(meta), path("${prefix}.metrics.txt"), emit: metrics
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    // Additional known UMIs may also be supplied inline with --umis via task.ext.args
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_corrected"
    def umis_command = umis ? "--umi-files ${umis}" : ''
    def rejects_command = keep_rejected ? "--rejects ${prefix}.rejects.bam" : ''
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    fgumi \\
        correct \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --min-distance ${min_distance} \\
        --metrics ${prefix}.metrics.txt \\
        --threads ${task.cpus} \\
        ${umis_command} \\
        ${rejects_command} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_corrected"
    def rejects_command = keep_rejected ? "touch ${prefix}.rejects.bam" : ''
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.metrics.txt
    ${rejects_command}
    """
}
