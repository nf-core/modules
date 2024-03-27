process CRABZ_DECOMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabz:0.9.0':
        'biocontainers/crabz:0.9.0' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("*.*"), emit: file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${archive}"
    """
    crabz \\
        $args \\
        -p $task.cpus \\
        -o ${prefix} \\
        $archive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabz: \$(crabz --version |& sed 's/[^:]*://')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${archive}"
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabz: \$(crabz --version |& sed 's/[^:]*://')
    END_VERSIONS
    """
}
