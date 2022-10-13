process LAST_LASTDB {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    def container_image = "/last:1250--h2e03b76_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("lastdb"), emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir lastdb
    lastdb \\
        $args \\
        -P $task.cpus \\
        lastdb/${prefix} \\
        $fastx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """
}
