process WISECONDORX_CONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.1--py39h67e14b5_0' :
        'biocontainers/wisecondorx:1.2.1--py39h67e14b5_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    WisecondorX \\
        convert \\
        ${bam} \\
        ${prefix}.npz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: \$(WisecondorX --version 2>&1 | sed 's/WisecondorX //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: \$(WisecondorX --version 2>&1 | sed 's/WisecondorX //')
    END_VERSIONS
    """
}