process WISECONDORX_CONVERT {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def VERSION = '1.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    WisecondorX convert \\
        ${bam} \\
        ${prefix}.npz \\
        ${reference} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """
}
