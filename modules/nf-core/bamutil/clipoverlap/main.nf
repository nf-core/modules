process BAMUTIL_CLIPOVERLAP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bamutil:1.0.15--h2e03b76_1'
        : 'biocontainers/bamutil:1.0.15--h2e03b76_1'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: stats_log, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_clipoverlap"
    """
    bam \\
        clipOverlap \\
        --in ${bam} \\
        --out ${prefix}.bam \\
        --stats \\
        ${args} \\
        2> >( tee ${prefix}.log >&2 ) \\
        | tee ${prefix}.err


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamutil: \$(bam clipOverlap 2>&1 | grep '^Version:' | sed 's/^Version: //;s/;.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_clipoverlap"
    """
    touch ${prefix}.bam
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamutil: \$(bam clipOverlap 2>&1 | grep '^Version:' | sed 's/^Version: //;s/;.*//')
    END_VERSIONS
    """
}
