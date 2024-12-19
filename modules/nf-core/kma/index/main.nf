process KMA_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.15--he4a0461_0' :
        'biocontainers/kma:1.4.15--he4a0461_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.kmaindex.*"),   emit: db
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}.kmaindex"
    """
    kma \\
        index \\
        -i ${fasta} \\
        -o ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma_index -v 2>&1) | sed 's/^KMA_index-\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}.kmaindex"
    """
    touch ${prefix}.comp.b
    touch ${prefix}.length.b
    touch ${prefix}.name
    touch ${prefix}.seq.b

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma_index -v 2>&1) | sed 's/^KMA_index-\$//')
    END_VERSIONS
    """
}
