process PIGZ_UNCOMPRESS {
    tag '$file'
    label 'process_low'

    conda "conda-forge::pigz"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE': //TODO add when multicontainer has been built
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    path zip

    output:
    path $zip.dropRight(3), emit: file
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    unpigz \\
        -p $task.cpus \\
        $args \\
        $zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch ${zip.dropRight(3)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\w*//' ))
    END_VERSIONS
    """
}
