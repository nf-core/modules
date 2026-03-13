process RATTLE_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ecoflowucl/rattle:v1.0' :
        'ecoflowucl/rattle:v1.0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("clusters.out"), emit: clusters
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def RATTLE_VERSION = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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
    def RATTLE_VERSION = "v1.0"
    """
    touch clusters.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rattle: $RATTLE_VERSION
    END_VERSIONS
    """
}
