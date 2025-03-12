process BBMAP_PILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bbmap_samtools_pigz:2a066f0214cc5eb0' :
        'community.wave.seqera.io/library/bbmap_samtools_pigz:79703e935236b43b' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats.txt"), emit: covstats
    tuple val(meta), path("*.hist.txt") , emit: hist
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pileup.sh \\
        -Xmx${task.memory.toGiga()}g \\
        in=${bam} \\
        out=${prefix}.coverage.stats.txt \\
        hist=${prefix}.coverage.hist.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
