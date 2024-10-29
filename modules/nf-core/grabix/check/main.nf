process GRABIX_CHECK {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grabix:0.1.8--hdcf5f25_9':
        'biocontainers/grabix:0.1.8--hdcf5f25_9' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), stdout, emit: compress_bgzip
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    grabix check ${input} | tr -d '\\n'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grabix: \$(grabix | sed -n -E 's/version: (.*)/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    \$(echo yes)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grabix: \$(grabix | sed -n -E 's/version: (.*)/\\1/p')
    END_VERSIONS
    """
}
