process GFATOOLS_STAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfatools:0.5--h577a1d6_5':
        'biocontainers/gfatools:0.5--h577a1d6_5' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('gfatools'), eval("gfatools version | sed '1!d; s/.* //'"), topic: versions, emit: versions_gfatools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gfatools \\
        stat \\
        $args \\
        $gfa \\
        > ${prefix}.stats
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    """
}
