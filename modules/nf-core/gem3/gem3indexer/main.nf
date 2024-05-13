process GEM3_GEM3INDEXER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gem3-mapper:3.6.1--h9d449c0_12':
        'biocontainers/gem3-mapper:3.6.1--h9d449c0_12' }"

    input:
    tuple val(meta), path(fasta)


    output:
    tuple val(meta), path("*.gem") , emit: index
    tuple val(meta), path("*.info"), emit: info
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gem-indexer \\
        -i ${fasta} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem3-indexer: \$(echo \$(gem-indexer --version 2>&1) | sed 's/v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.gem
    touch ${prefix}.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem3-indexer: \$(echo \$(gem-indexer --version 2>&1) | sed 's/v//')
    END_VERSIONS
    """
}
