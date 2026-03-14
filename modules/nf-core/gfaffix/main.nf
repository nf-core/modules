process GFAFFIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfaffix:0.2.1--hc1c3326_0' :
        'biocontainers/gfaffix:0.2.1--hc1c3326_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    tuple val(meta), path("*.txt"), emit: affixes
    tuple val("${task.process}"), val('gfaffix'), eval('gfaffix --version'), emit: versions_gfaffix, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gfaffix \\
        $args \\
        $gfa \\
        -o ${prefix}.gfaffix.gfa > ${prefix}.affixes.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gfaffix.gfa
    touch ${prefix}.affixes.txt
    """
}
