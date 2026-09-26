process SENTIEON_COLLECTVCMETRICS {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7ee64f3b4cd58eaa6ed3a7a0769c0a5f4fbd150842fa64358831970bacaf37e6/data'
        : 'community.wave.seqera.io/library/sentieon:202503.03--df1987151f8b6d33'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(dbsnp), path(dbsnp_tbi)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), path(interval)

    output:
    tuple val(meta), path("*.variant_calling_detail_metrics"),  emit: metrics
    tuple val(meta), path("*.variant_calling_summary_metrics"), emit: summary
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def interval_cmd = interval ? "--interval ${interval}" : ""
    """
    sentieon \\
        driver \\
        -t ${task.cpus} \\
        -r ${fasta} \\
        ${interval_cmd} \\
        ${args} \\
        --algo CollectVCMetrics \\
        -v ${vcf} \\
        -d ${dbsnp} \\
        ${args2} \\
        ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_calling_detail_metrics
    touch ${prefix}.variant_calling_summary_metrics
    """
}
