process PIGZ_UNCOMPRESS {
    label 'process_low'
    //stageInMode 'copy' // this directive can be set in case the original input should be kept

    conda "conda-forge::pigz"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8':
        'biocontainers/pigz:2.8' }"

    input:
    path zip

    output:
    path "${uncompressed_filename}" , emit: file
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    uncompressed_filename = zip.toString() - '.gz'
    // calling pigz -f to make it follow symlinks
    """
    unpigz \\
        -p $task.cpus \\
        -fk \\
        $args \\
        ${zip}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    uncompressed_filename = zip.toString() - '.gz'
    """
    touch ${zip.dropRight(3)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\w*//' ))
    END_VERSIONS
    """
}
