process ANGSD_DOCOUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.939--h468462d_0':
        'biocontainers/angsd:0.939--h468462d_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(minqfile)

    output:
    tuple val(meta), path("*.depthSample"), emit: depth_sample, optional: true
    tuple val(meta), path("*.depthGlobal"), emit: depth_global, optional: true
    tuple val(meta), path("*.qs")         , emit: qs          , optional: true
    tuple val(meta), path("*.pos.gz")     , emit: pos         , optional: true
    tuple val(meta), path("*.counts.gz")  , emit: counts      , optional: true
    tuple val(meta), path("*.icnts.gz")   , emit: icounts     , optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def minq   = minqfile ? "-minQfile ${minqfile}" : ""
    """
    ls -1 *.bam > bamlist.txt

    angsd \\
        -nThreads ${task.cpus} \\
        -doCounts 1 \\
        ${args} \\
        -bam bamlist.txt \\
        -out ${prefix} \\
        ${minq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    depthSample_cmd = args.contains('-doDepth')       ? "touch ${prefix}.depthSample"       : ''
    depthGlobal_cmd = args.contains('-doDepth')       ? "touch ${prefix}.depthGlobal"       : ''
    qs_cmd          = args.contains('-doQsDist')      ? "touch ${prefix}.qs"                : ''
    pos_cmd         = args.contains('-dumpCounts')    ? "echo | gzip > ${prefix}.pos.gz"    : ''
    counts_cmd      = args =~ /-dumpCounts\s+[12345]/ ? "echo | gzip > ${prefix}.counts.gz" : ''
    icounts_cmd     = args.contains('-iCounts')       ? "echo | gzip > ${prefix}.icnts.gz"  : ''

    """
    ${depthSample_cmd}
    ${depthGlobal_cmd}
    ${qs_cmd}
    ${pos_cmd}
    ${counts_cmd}
    ${icounts_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}
