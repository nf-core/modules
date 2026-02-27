process COOLER_ZOOMIFY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*.mcool"), emit: mcool
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions_cooler, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mcool
    """
}
