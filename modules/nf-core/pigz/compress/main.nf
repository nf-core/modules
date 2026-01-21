process PIGZ_COMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8':
        'biocontainers/pigz:2.8' }"

    input:
    tuple val(meta), path(raw_file)

    output:
    tuple val(meta), path("$archive"), emit: archive
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    archive = raw_file.toString() + ".gz"
    """
    # Note: needs --stdout for pigz to avoid the following issue:
    #   pigz: skipping: ${raw_file} is a symbolic link
    pigz --processes $task.cpus --stdout --force ${args} ${raw_file} > ${archive}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz[[:space:]]*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    archive = raw_file.toString() + ".gz"
    """
    touch ${archive}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz[[:space:]]*//' )
    END_VERSIONS
    """
}
