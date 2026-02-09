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
    tuple val("${task.process}"), val('xz'), eval('xz --version | sed \'1!d;s/\\([^0-9.]*\\)//g\''), topic: versions, emit: versions_xz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    archive = raw_file.toString() + ".xz"
    """
    # needs --stdout for xz to avoid the following issue:
    # xz: ${raw_file}: Is a symbolic link, skipping
    xz -T $task.cpus --stdout ${args} ${raw_file} > ${archive}
    """

    stub:
    def args = task.ext.args ?: ''
    archive = raw_file.toString() + ".xz"
    """
    echo "${args}"
    touch "${archive}"
    """
}
