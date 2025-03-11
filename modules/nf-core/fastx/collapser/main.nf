process FASTX_COLLAPSER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastx_toolkit:0.0.14--hdbdd923_11':
        'biocontainers/fastx_toolkit:0.0.14--hdbdd923_11' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastx_collapser \\
        $args \\
        -i $fastx \\
        -o ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastx:  \$(echo \$(fastx_collapser -h) | sed -nE 's/.*([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastx:  \$(echo \$(fastx_collapser -h) | sed -nE 's/.*([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' ))
    END_VERSIONS
    """
}
