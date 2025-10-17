process SPRING_COMPRESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spring:1.1.1--h9f5acd7_0' :
        'biocontainers/spring:1.1.1--h9f5acd7_0' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta), path("*.spring"), emit: spring
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def input = meta.single_end ? "-i ${fastq1}" : "-i ${fastq1} ${fastq2}"

    """
    spring \\
        -c \\
        -g \\
        -t ${task.cpus} \\
        $args \\
        ${input} \\
        -o ${prefix}.spring

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spring: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.spring

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spring: ${VERSION}
    END_VERSIONS
    """

}
