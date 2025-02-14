process SPLITUBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tagbam:0.1.0--hc9368f3_0':
        'biocontainers/tagbam:0.1.0-hc9368f3_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    tagbam \\
        $args \\
        --threads $task.cpus \\
        --input $bam \\
        --output-file {prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splitubam: \$(splitubam --version | sed 's/splitubam //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splitubam: \$(splitubam --version | sed 's/splitubam //')
    END_VERSIONS
    """
}
