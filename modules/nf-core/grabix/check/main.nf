process GRABIX_CHECK {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grabix:0.1.8--hdcf5f25_9':
        'biocontainers/grabix:0.1.8--hdcf5f25_9' }"

    input:
    path(input)

    output:
    env COMPRESS_BGZIP  , emit: compress_bgzip
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    COMPRESS_BGZIP=\$(grabix check ${input})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grabix: \$(echo \$(grabix | sed -n -E 's/version: (.*)/\1/p'))
    END_VERSIONS
    """
}
