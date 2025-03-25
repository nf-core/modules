process POLYPOLISH {
    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/polypolish:0.6.0--hdbdd923_0':
        'biocontainers/polypolish:0.6.0' }"

    input:
    // input:  polypolish polish [OPTIONS] <ASSEMBLY> [SAM]
    tuple val(meta), path(fasta)
    tuple val(meta), path(sam)

    output:
    // polypolish polish draft.fasta filtered_1.sam filtered_2.sam > polished.fasta
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    polypolish \\
        polish \\
        $args \\
        ${fasta} \\
        $sam > ${prefix}_polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(polypolish polish --version |& sed '1!d ; s/polypolish polish //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(polypolish polish --version |& sed '1!d ; s/polypolish polish //')
    END_VERSIONS
    """
}
