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
    tuple val(meta), path("trim/$prefix/*.gfa"),  emit: gfa
    tuple val(meta), path("trim/$prefix/*.yaml"), emit: stats
    tuple val("${task.process}"), val("autocycler"), eval("autocycler --version |  sed 's/^[^ ]* //'"), emit: versions_autocycler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    autocycler trim \\
        $args \\
        --threads $task.cpus \\
        -c .

    mkdir -p trim/$prefix
    mv 2_trimmed.{gfa,yaml} trim/$prefix
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """

    mkdir -p trim/$prefix
    touch trim/$prefix/2_trimmed.gfa
    touch trim/$prefix/2_trimmed.yaml
    """
}
