process RCLONE_COPY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5dfd28fd0090c69f57c9bd93ea3235d8df194e8f269b5cb3027b6b59bff567d5/data'
            : 'community.wave.seqera.io/library/rclone:1.74.3--2ef33c5b9132aa97' }"

    input:
    tuple val(meta), val(source_path), val(destination_path), path(filter_file)
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
        def matcher = (source_string =~ /^(https?:\/\/[^\/]+)(\/.*)?$/)
        if (!matcher.matches()) {
            throw new IllegalArgumentException("Invalid HTTP(S) source '${source_string}' for sample '${meta.id}'.")
        }
        http_url_arg = "--http-url '${matcher[0][1]}'"
        rclone_source = ":http:${(matcher[0][2] ?: '/').replaceFirst('^/', '')}"
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
