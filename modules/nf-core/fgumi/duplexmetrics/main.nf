process FGUMI_DUPLEXMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a510706f3481fae12ff6100d6e4ad298b8bf464a2d93a6afe35e9cf26542d080/data':
        'community.wave.seqera.io/library/fgumi:0.2.0--fe028e7a64e5da27' }"

    input:
    tuple val(meta), path(bam), path(intervals)

    output:
    tuple val(meta), path("*.family_sizes.txt"), emit: family_sizes
    tuple val(meta), path("*.duplex_family_sizes.txt"), emit: duplex_family_sizes
    tuple val(meta), path("*.duplex_yield_metrics.txt"), emit: duplex_yield_metrics
    tuple val(meta), path("*.umi_counts.txt"), emit: umi_counts
    tuple val(meta), path("*.duplex_qc.pdf"), emit: duplex_qc, optional: true
    tuple val(meta), path("*.duplex_umi_counts.txt"), emit: duplex_umi_counts, optional: true
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals ${intervals}" : ''
    """
    fgumi \\
        duplex-metrics \\
        --input ${bam} \\
        --output ${prefix} \\
        ${intervals_command} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_duplex_umi = args.contains("--duplex-umi-counts") ? "touch ${prefix}.duplex_umi_counts.txt" : ""
    """
    touch ${prefix}.family_sizes.txt
    touch ${prefix}.duplex_family_sizes.txt
    touch ${prefix}.duplex_yield_metrics.txt
    touch ${prefix}.umi_counts.txt
    touch ${prefix}.duplex_qc.pdf
    ${touch_duplex_umi}
    """
}
