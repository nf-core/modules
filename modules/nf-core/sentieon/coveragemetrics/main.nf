process SENTIEON_COVERAGEMETRICS {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a64461f38d76bebea8e21441079e76e663e1168b0c59dafee6ee58440ad8c8ac/data' :
        'community.wave.seqera.io/library/sentieon:202308.03--59589f002351c221' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(interval)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), path(gene_list)

    output:
    tuple val(meta), path("$prefix")                                                       , optional: true, emit: per_locus
    tuple val(meta), path("${prefix}.${partitions_output}_summary")                        , optional: true, emit: sample_summary
    tuple val(meta), path("${prefix}.${partitions_output}_interval_statistics")            , optional: true, emit: statistics
    tuple val(meta), path("${prefix}.${partitions_output}_cumulative_coverage_counts")     , optional: true, emit: coverage_counts
    tuple val(meta), path("${prefix}.${partitions_output}_cumulative_coverage_proportions"), optional: true, emit: coverage_proportions
    tuple val(meta), path("${prefix}.${partitions_output}_interval_summary")               , optional: true, emit: interval_summary
    path "versions.yml"                                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix  = task.ext.prefix ?: "${meta.id}"
    def args           = task.ext.args  ?: ''
    def args2          = task.ext.args2 ?: ''
    def input          = bam.sort().collect{"-i $it"}.join(' ')
    def interval_cmd   = interval   ? "--interval $interval"                : ""
    def gene_list_cmd  = gene_list  ? "--gene_list ${gene_list}"            : ""
    // Glob that matches any version of 'sample_library_platform_center'.
    partitions_output = "{sample,}{_library,}{_platform,}{_center,}{_readgroup,}"
    """
    sentieon \\
        driver \\
        -t $task.cpus \\
        -r $fasta \\
        $interval_cmd \\
        $input \\
        $args \\
        --algo CoverageMetrics \\
        $args2 \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}
    touch ${prefix}.sample_interval_statistics
    touch ${prefix}.sample_cumulative_coverage_counts
    touch ${prefix}.sample_cumulative_coverage_proportions
    touch ${prefix}.sample_interval_summary
    touch ${prefix}.sample_cumulative_coverage_counts
    touch ${prefix}.sample_cumulative_coverage_proportions
    touch ${prefix}.sample_interval_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
