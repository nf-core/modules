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
    tuple val(meta), path("*rclone-copy.log"), emit: log
    tuple val("${task.process}"), val('rclone'), eval("rclone --version | sed -n '1s/^rclone v//p'"), topic: versions, emit: versions_rclone

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def configArg = rclone_config ? "--config '${rclone_config}'" : ''
    def transfers = Math.max(1, task.cpus.intdiv(2))
    def checkers = task.cpus

    // Handle HTTP URLs: split into --http-url base and :http:relative_path
    def source_string = source_path.toString()
    def rclone_source
    def http_url_arg = ''

    if (source_string ==~ /^https?:\/\/.*/) {
        def matcher = (source_string =~ /^(https?:\/\/[^\/]+)(\/.*)$/)
        http_url_arg = "--http-url '${matcher[0][1]}'"
        rclone_source = ":http:${matcher[0][2].replaceFirst('^/', '')}"
    } else {
        rclone_source = source_string.replaceFirst('^([a-zA-Z][a-zA-Z0-9+.-]*)://', '$1:')
    }

    """
    rclone ${configArg} copy \\
        ${http_url_arg} \\
        ${args} \\
        --log-file "${meta.id}-rclone-copy.log" \\
        --transfers ${transfers} \\
        --checkers ${checkers} \\
        "${rclone_source}" \\
        "${destination_path}"
    """

    stub:
    """
    touch rclone-copy.log
    """
}
