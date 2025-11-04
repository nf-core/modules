process AUTOCYCLER_TRIM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler:0.5.2--h3ab6199_0':
        'biocontainers/autocycler:0.5.2--h3ab6199_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.gfa"),  emit: gfa
    tuple val(meta), path("*.yaml"), emit: stats
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    """
    autocycler trim \\
        $args \\
        --threads $task.cpus \\
        -c .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    """

    touch 2_trimmed.gfa
    touch 2_trimmed.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """
}
