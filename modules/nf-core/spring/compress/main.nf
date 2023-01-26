process SPRING_COMPRESS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::spring=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spring:1.1.0--h9f5acd7_0' :
        'quay.io/biocontainers/spring:1.1.0--h9f5acd7_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.spring"), emit: spring
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    if (meta.single_end) {
        """
        spring -c -g -t ${task.cpus} $args -i '${fastq}' -o '${prefix}.spring'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            spring: ${VERSION}
        END_VERSIONS
        """
    } else {
        """
        spring -c -g -t ${task.cpus} $args -i '${fastq[0]}' '${fastq[1]}' -o '${prefix}.spring'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            spring: ${VERSION}
        END_VERSIONS
        """
    }
}
