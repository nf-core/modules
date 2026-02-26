process PURGEDUPS_PURGEDUPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--py39h7132678_1':
        'biocontainers/purge_dups:1.2.6--py39h7132678_1' }"

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
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version
    """
    purge_dups \\
        $args \\
        -T $cutoff \\
        -c $basecov \\
        $paf > ${prefix}.dups.bed 2>| >(tee ${prefix}.purge_dups.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version
    """
    touch ${prefix}.dups.bed
    touch ${prefix}.purge_dups.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: ${VERSION}
    END_VERSIONS
    """
}
