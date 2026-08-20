process FGUMI_SIMPLEXMETRICS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e6/e613097ca7c84595a8683b6b042d113ac40e1936a809477aa95fcd1f8f3bfca2/data'
:         'community.wave.seqera.io/library/fgumi:0.6.0--c97194d17da0d1cd' }"

    input:
    tuple val(meta), path(bam), path(intervals)

    output:
    tuple val(meta), path("*.family_sizes.txt"), emit: family_sizes
    tuple val(meta), path("*.simplex_yield_metrics.txt"), emit: simplex_yield_metrics
    tuple val(meta), path("*.umi_counts.txt"), emit: umi_counts
    tuple val(meta), path("*.simplex_qc.pdf"), emit: pdf, optional: true
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals ${intervals}" : ''
    """
    fgumi \\
        simplex-metrics \\
        --input ${bam} \\
        --output ${prefix} \\
        ${intervals_command} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}.family_sizes.txt
    touch ${prefix}.simplex_yield_metrics.txt
    touch ${prefix}.umi_counts.txt
    touch ${prefix}.simplex_qc.pdf
    """
}
