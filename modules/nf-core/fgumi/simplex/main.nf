process FGUMI_SIMPLEX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a510706f3481fae12ff6100d6e4ad298b8bf464a2d93a6afe35e9cf26542d080/data'
        : 'community.wave.seqera.io/library/fgumi:0.2.0--fe028e7a64e5da27'}"

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
    prefix = task.ext.prefix ?: "${meta.id}_simplex_unmapped"
    def rejects_command = keep_rejected ? "--rejects ${prefix}.rejects.bam" : ''

    if ("${grouped_bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    fgumi simplex \\
        --input ${grouped_bam} \\
        --output ${prefix}.bam \\
        --min-reads ${min_reads} \\
        --threads ${task.cpus} \\
        --stats ${prefix}.stats.txt \\
        ${rejects_command} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_simplex_unmapped"
    if ("${grouped_bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.rejects.bam
    touch ${prefix}.stats.txt
    """
}
