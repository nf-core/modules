process PURGEDUPS_SPLITFA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--py39h7132678_1':
        'biocontainers/purge_dups:1.2.6--py39h7132678_1' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.split.fasta.gz"), emit: split_fasta
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version
    """
    split_fa $args $assembly | gzip -c > ${prefix}.split.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version
    """
    echo | gzip > ${prefix}.split.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: $VERSION
    END_VERSIONS
    """

}
