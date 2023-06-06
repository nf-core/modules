process SHASUM {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

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
