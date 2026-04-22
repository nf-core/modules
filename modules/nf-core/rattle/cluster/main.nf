process RATTLE_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rattle:1.0--h5ca1c30_0' :
        'biocontainers/rattle:1.0--h5ca1c30_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("clusters.out"), emit: clusters
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def RATTLE_VERSION = "1.0"
    """
    rattle \\
        cluster $args \\
        -t $task.cpus \\
        -i $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rattle: $RATTLE_VERSION
    END_VERSIONS
    """

    stub:
    def RATTLE_VERSION = "1.0"
    """
    touch clusters.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rattle: $RATTLE_VERSION
    END_VERSIONS
    """
}
