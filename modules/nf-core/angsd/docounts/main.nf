process ANGSD_DOCOUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::angsd=0.939"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.939--h468462d_0':
        'biocontainers/angsd:0.939--h468462d_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(minqfile)

    output:
    tuple val(meta), path("*.depthSample"), optional: true, emit: depth_sample
    tuple val(meta), path("*.depthGlobal"), optional: true, emit: depth_global
    tuple val(meta), path("*.qs")         , optional: true, emit: qs
    tuple val(meta), path("*.pos.gz")     , optional: true, emit: pos
    tuple val(meta), path("*.counts.gz")  , optional: true, emit: counts
    tuple val(meta), path("*.icnts.gz")   , optional: true, emit: icounts
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def minq = minqfile ? "-minQfile ${minqfile}" : ""
    """
    ls -1 *.bam > bamlist.txt

    angsd \\
        -nThreads ${task.cpus} \\
        -doCounts 1 \\
        $args \\
        -bam bamlist.txt \\
        -out ${prefix} \\
        $minq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}
