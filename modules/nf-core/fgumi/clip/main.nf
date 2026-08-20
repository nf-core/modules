process FGUMI_CLIP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e6/e613097ca7c84595a8683b6b042d113ac40e1936a809477aa95fcd1f8f3bfca2/data':
        'community.wave.seqera.io/library/fgumi:0.6.0--c97194d17da0d1cd' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics, optional: true
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_clipped"
    // fgumi clip cannot combine --metrics with --threads, so only request threads when metrics are not
    def threads = args.contains('--metrics') ? '' : "--threads ${task.cpus}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    fgumi \\
        clip \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --reference ${fasta} \\
        ${threads} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_clipped"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    """
}
