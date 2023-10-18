process PIGZ_COMPRESS {
    tag '$bam'
    label 'process_low'

    conda "conda-forge::pigz"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE': //TODO add when multicontainer has been built
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    path file

    output:
    path "*.gz"        , emit: zip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    pigz \\
        -p $task.cpus \\
        $args \\
        $file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch ${file}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\w*//' ))
    END_VERSIONS
    """
}
