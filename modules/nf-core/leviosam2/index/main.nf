process LEVIOSAM2_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leviosam2:0.4.2--h4ac6f70_0':
        'biocontainers/leviosam2:0.4.2--h4ac6f70_0' }"

    input:
    tuple val(meta), path(fai)
    path(chain)

    output:
    tuple val(meta), path("*.clft"), emit: clft
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    leviosam2 \\
        index \\
        -c ${chain} \\
        -p ${prefix} \\
        -F ${fai}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leviosam2: \$(leviosam2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clft

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leviosam2: \$(leviosam2 --version)
    END_VERSIONS
    """
}
