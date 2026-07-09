process RCLONE_COPY {

    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c947c1a7171daf074310295417d0f0afe879275e1543fa5fc2a9711e7c2c72ab/data'
        : 'community.wave.seqera.io/library/rclone:1.65.0--ff88b2e0040147be'}"

    input:
    tuple val(meta), val(source_path), val(destination_path)
    path rclone_config

    output:
    tuple val(meta), path("rclone-copy.log"), emit: log
    tuple val("${task.process}"), val('rclone'), eval("rclone version | head -n1 | sed 's/rclone v//'"), topic: versions, emit: versions_rclone

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def configArg = rclone_config ? "--config '${rclone_config}'" : ''
    def transfers = task.cpus
    def checkers = task.cpus

    """
    rclone ${configArg} copy ${args} \\
        --log-file rclone-copy.log \\
        --transfers ${transfers} \\
        --checkers ${checkers} \\
        "${source_path}" \\
        "${destination_path}"
    """

    stub:
    """
    touch rclone-copy.log
    """
}
