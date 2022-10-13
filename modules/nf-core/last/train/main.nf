process LAST_TRAIN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    def container_image = "/last:1250--h2e03b76_0"
                                             container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fastx)
    path  index

    output:
    tuple val(meta), path("*.par"), emit: param_file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)

    last-train \\
        $args \\
        -P $task.cpus \\
        ${index}/\$INDEX_NAME \\
        $fastx \\
        > ${prefix}.\$INDEX_NAME.par

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version | sed 's/lastdb //')
    END_VERSIONS
    """
}
