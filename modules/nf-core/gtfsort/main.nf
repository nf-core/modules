process GTFSORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gtfsort:0.2.2--h4ac6f70_0'
        : 'biocontainers/gtfsort:0.2.2--h4ac6f70_0'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.sorted.${gtf.extension}"), emit: gtf
    tuple val("${task.process}"), val('gtfsort'), eval('gtfsort --version |& sed "s/gtfsort //"'), emit: versions_gtfsort, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gtfsort \\
        -i ${gtf} \\
        -o ${prefix}.sorted.${gtf.extension} \\
        ${args} \\
        -t ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.sorted.${gtf.extension}
    """
}
