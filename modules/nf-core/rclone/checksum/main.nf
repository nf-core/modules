process RCLONE_CHECKSUM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5dfd28fd0090c69f57c9bd93ea3235d8df194e8f269b5cb3027b6b59bff567d5/data'
            : 'community.wave.seqera.io/library/rclone:1.74.3--2ef33c5b9132aa97' }"

    input:
    tuple val(meta), path(sumfile), val(hash), val(destination)
    path rclone_config

    output:
    tuple val(meta), path("${prefix}.combined.txt")       , emit: combined      , optional: true
    tuple val(meta), path("${prefix}.differ.txt")         , emit: differ        , optional: true
    tuple val(meta), path("${prefix}.missing_on_dst.txt") , emit: missing_on_dst, optional: true
    tuple val(meta), path("${prefix}.missing_on_src.txt") , emit: missing_on_src, optional: true
    tuple val(meta), path("${prefix}.match.txt")          , emit: match         , optional: true
    tuple val(meta), path("${prefix}.error.txt")          , emit: error         , optional: true
    tuple val("${task.process}"), val('rclone'), eval("rclone --version | sed -n '1s/^rclone v//p'"), topic: versions, emit: versions_rclone

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def configArg = rclone_config ? "--config ${rclone_config}" : ''

    """
    rclone checksum ${configArg} \\
        --copy-links \\
        $args \\
        --combined ${prefix}.combined.txt \\
        --differ ${prefix}.differ.txt \\
        --missing-on-dst ${prefix}.missing_on_dst.txt \\
        --missing-on-src ${prefix}.missing_on_src.txt \\
        --match ${prefix}.match.txt \\
        --error ${prefix}.error.txt \\
        --checkers $task.cpus \\
        $hash \\
        $sumfile \\
        ${destination} || true

    # Do not emit empty output files
    for f in *.txt; do
        [ -s "\$f" ] || rm -f "\$f"
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch \\
        ${prefix}.combined.txt \\
        ${prefix}.differ.txt \\
        ${prefix}.missing_on_dst.txt \\
        ${prefix}.missing_on_src.txt \\
        ${prefix}.match.txt \\
        ${prefix}.error.txt
    """
}
