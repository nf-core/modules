process SENTIEON_COVERAGEMETRICS {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

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
    touch ${prefix}.sample_interval_statistics
    touch ${prefix}.sample_cumulative_coverage_counts
    touch ${prefix}.sample_cumulative_coverage_proportions
    touch ${prefix}.sample_interval_summary
    touch ${prefix}.sample_cumulative_coverage_counts
    touch ${prefix}.sample_cumulative_coverage_proportions
    touch ${prefix}.sample_interval_summary
    """
}
