
process MIRTOP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/mirtop_pybedtools_pysam_samtools_pruned:60b8208f3dbb2910"

    input:
    tuple val(meta), path(mirtop_gff)

    output:
    tuple val(meta), path("stats/*.txt")        , emit: txt
    tuple val(meta), path("stats/*_stats.log")  , emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mirtop \\
        stats \\
        $args \\
        --out stats \\
        $mirtop_gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir stats
    touch stats/${prefix}.txt
    touch stats/${prefix}_stats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """
}
