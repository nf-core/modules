process FGUMI_DUPLEX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64e8f594b6f0dd879bc5abbe4ca70b6b761e1920e407d9e1c7d27b89004aac34/data'
        : 'community.wave.seqera.io/library/fgumi:0.5.0--a2d14bf52f73eaef'}"

    input:
    tuple val(meta), path(grouped_bam)
    val min_reads
    val keep_rejected

    output:
    tuple val(meta), path("${prefix}.bam")        , emit: bam
    tuple val(meta), path("${prefix}.rejects.bam"), emit: rejects, optional: true
    tuple val(meta), path("${prefix}.stats.txt")  , emit: stats
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_duplex_unmapped"
    def rejects_command = keep_rejected ? "--rejects ${prefix}.rejects.bam" : ''

    if ("${grouped_bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    fgumi duplex \\
        --input ${grouped_bam} \\
        --output ${prefix}.bam \\
        --min-reads ${min_reads} \\
        --threads ${task.cpus} \\
        --stats ${prefix}.stats.txt \\
        ${rejects_command} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_duplex_unmapped"
    def rejects_command = keep_rejected ? "touch ${prefix}.rejects.bam" : ''
    if ("${grouped_bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    ${rejects_command}
    touch ${prefix}.stats.txt
    """
}
