process FLASH {
    tag "$meta.id"
    label 'process_medium'

        def container_image = "/flash:1.2.11--hed695b0_5"
                                        container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    flash \\
        $args \\
        -o ${prefix} \\
        -z \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flash: \$(echo \$(flash --version 2>&1) | sed 's/^.*FLASH v//; s/ .*\$//')
    END_VERSIONS
    """
}
