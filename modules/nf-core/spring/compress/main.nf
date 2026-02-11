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
    tuple val("${task.process}"), val('spring'), val("$VERSION"), topic: versions, emit: versions_spring

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '1.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def input = meta.single_end ? "-i ${fastq1}" : "-i ${fastq1} ${fastq2}"

    """
    spring \\
        -c \\
        -g \\
        -t ${task.cpus} \\
        $args \\
        ${input} \\
        -o ${prefix}.spring
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '1.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.spring
    """

}
