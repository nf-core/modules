process CIRCULARMAPPER_CIRCULARGENERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circularmapper:1.93.5--h4a94de4_1':
        'biocontainers/circularmapper:1.93.5--h4a94de4_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_${task.ext.prefix}.${fasta.extension}") , emit: fasta
    tuple val(meta), path("*_${task.ext.prefix}_elongated")          , emit: elongated
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.93.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    circulargenerator \\
        -Xmx${task.memory.toGiga()}g \\
        $args \\
        -i $fasta
    mv ${fasta}_${prefix}_elongated ${fasta.baseName}_${prefix}.${fasta.extension}_${prefix}_elongated

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CircularMapper: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.93.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${fasta.baseName}_${prefix}.${fasta.extension}
    touch ${fasta.baseName}_${prefix}.${fasta.extension}_${prefix}_elongated

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CircularMapper: ${VERSION}
    END_VERSIONS
    """
}
