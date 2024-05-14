process COOLER_ZOOMIFY {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::cooler=0.9.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.9.2--pyh7cba7a3_0' :
        'biocontainers/cooler:0.9.2--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*.mcool"), emit: mcool
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooler zoomify \\
        $args \\
        -n $task.cpus \\
        -o ${prefix}.mcool \\
        $cool

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mcool

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
