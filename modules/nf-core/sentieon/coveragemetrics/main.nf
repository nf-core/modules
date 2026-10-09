process SENTIEON_COVERAGEMETRICS {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7ee64f3b4cd58eaa6ed3a7a0769c0a5f4fbd150842fa64358831970bacaf37e6/data'
        : 'community.wave.seqera.io/library/sentieon:202503.03--df1987151f8b6d33'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(interval)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), path(gene_list)

    output:
    tuple val(meta), path("${prefix}"),                                                      emit: per_locus,            optional: true
    tuple val(meta), path("${prefix}.${partitions_output}_summary"),                         emit: sample_summary,       optional: true
    tuple val(meta), path("${prefix}.${partitions_output}_interval_statistics"),             emit: statistics,           optional: true
    tuple val(meta), path("${prefix}.${partitions_output}_cumulative_coverage_counts"),      emit: coverage_counts,      optional: true
    tuple val(meta), path("${prefix}.${partitions_output}_cumulative_coverage_proportions"), emit: coverage_proportions, optional: true
    tuple val(meta), path("${prefix}.${partitions_output}_interval_summary"),                emit: interval_summary,     optional: true
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def input = bam.sort().collect {in -> "-i ${in}" }.join(' ')
    def interval_cmd = interval ? "--interval ${interval}" : ""
    def gene_list_cmd = gene_list ? "--gene_list ${gene_list}" : ""
    // Glob that matches any version of 'sample_library_platform_center'.
    partitions_output = "{sample,}{_library,}{_platform,}{_center,}{_readgroup,}"
    """
    sentieon \\
        driver \\
        -t ${task.cpus} \\
        -r ${fasta} \\
        ${interval_cmd} \\
        ${gene_list_cmd} \\
        ${input} \\
        ${args} \\
        --algo CoverageMetrics \\
        ${args2} \\
        ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Glob that matches any version of 'sample_library_platform_center'.
    partitions_output = "{sample,}{_library,}{_platform,}{_center,}{_readgroup,}"
    """
    touch ${prefix}
    touch ${prefix}.sample_summary
    touch ${prefix}.sample_interval_statistics
    touch ${prefix}.sample_cumulative_coverage_counts
    touch ${prefix}.sample_cumulative_coverage_proportions
    touch ${prefix}.sample_interval_summary
    """
}
