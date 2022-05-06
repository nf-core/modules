process MD5SUM {
    tag "$meta.id"
    label 'process_low'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using bcl-convert. Please use docker or singularity containers."
    }
    container "debian:bullseye-slim"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.md5"), emit: checksum
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    md5sum \\
        $args \\
        ${file} \\
        > ${file}.md5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        md5sum: \$(echo \$(md5sum --version 2>&1 | head -n 1| sed 's/^.*) //;' ))
    END_VERSIONS
    """
}
