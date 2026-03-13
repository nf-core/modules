process BBMAP_PILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats.txt"), emit: covstats
    tuple val(meta), path("*.hist.txt") , emit: hist
    tuple val("${task.process}"), val('bbmap'), eval("bbversion.sh | grep -v 'Duplicate cpuset'"), emit: versions_bbmap, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

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
    """
}
