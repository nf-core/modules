process SHASUM {
    tag "$meta.id"
    label 'process_low'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using bcl-convert. Please use docker or singularity containers."
    }
    container "debian:bullseye-slim"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.sha256"), emit: checksum
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sha256sum \\
        $args \\
        ${file} \\
        > ${file}.sha256

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sha256sum: \$(echo \$(sha256sum --version 2>&1 | head -n 1| sed 's/^.*) //;' ))
    END_VERSIONS
    """
}
