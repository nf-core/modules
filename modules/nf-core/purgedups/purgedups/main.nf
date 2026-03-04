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
    tuple val(meta), path("*.dups.bed.gz")   , emit: bed
    tuple val(meta), path("*.purge_dups.log"), emit: log
    // WARN: Incorrect version printed inside the container, please check this if bumping version ( \$( purge_dups -h |& sed '3!d; s/.*: //' ))
    tuple val("${task.process}"), val('purge_dups'), val('1.2.6'), emit: versions_purgedups, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    purge_dups \\
        ${args} \\
        -T ${cutoff} \\
        -c ${basecov} \\
        ${paf} | gzip \\
        > ${prefix}.dups.bed.gz \\
        2>| >(tee ${prefix}.purge_dups.log >&2)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    echo "" | gzip > ${prefix}.dups.bed.gz
    touch ${prefix}.purge_dups.log
    """
}
