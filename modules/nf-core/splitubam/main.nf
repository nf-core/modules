process SPLITUBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/splitubam:0.1.1--hc9368f3_0':
        'biocontainers/splitubam:0.1.1--hc9368f3_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    splitubam \\
        $args \\
        --threads $task.cpus \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splitubam: \$(splitubam --version | sed 's/splitubam //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def match = (args =~ /--split\s+(\d+)/)
    def create_cmd = ""
    if (match) {
        def n_splits = match[0][1].toInteger()
        (1..n_splits).each { i ->
            def formattedIteration = String.format('%03d', i)
            create_cmd += "touch ${formattedIteration}.${bam}.bam\n"
        }
    } else { error("No `--split N` detected in args") }
    """
    $create_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splitubam: \$(splitubam --version | sed 's/splitubam //')
    END_VERSIONS
    """
}
