process PURGEDUPS_PURGEDUPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8a/8afcc1222a0512a4fae7a7e7315ccfc841f3578df600b99d3dc563c3a8361352/data':
        'community.wave.seqera.io/library/purge_dups:1.2.6--1966ab26985f9f67' }"

    input:
    tuple val(meta), path(basecov), path(cutoff), path(paf)

    output:
    tuple val(meta), path("*.dups.bed")      , emit: bed
    tuple val(meta), path("*.purge_dups.log"), emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    purge_dups \\
        $args \\
        -T $cutoff \\
        -c $basecov \\
        $paf > ${prefix}.dups.bed 2> >(tee ${prefix}.purge_dups.log >&2)


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: \$( purge_dups -h |& sed '3!d; s/.*: //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.dups.bed
    touch ${prefix}.purge_dups.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: \$( purge_dups -h |& sed '3!d; s/.*: //' )
    END_VERSIONS
    """
}
