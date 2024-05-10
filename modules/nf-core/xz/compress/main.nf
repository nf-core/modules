process XZ_COMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-796b0610595ad1995b121d0b85375902097b78d4:a3a3220eb9ee55710d743438b2ab9092867c98c6-0' :
        'biocontainers/mulled-v2-796b0610595ad1995b121d0b85375902097b78d4:a3a3220eb9ee55710d743438b2ab9092867c98c6-0' }"

    input:
    tuple val(meta), path(raw_file)

    output:
    tuple val(meta), path("$archive"), emit: archive
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    archive = raw_file.toString() + ".xz"
    """
    # needs --stdout for xz to avoid the following issue:
    # xz: ${raw_file}: Is a symbolic link, skipping
    xz -T $task.cpus --stdout ${args} ${raw_file} > ${archive}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xz: \$(xz --version | head -n1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    archive = raw_file.toString() + ".xz"
    """
    touch "${archive}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xz: \$(xz --version | head -n1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
