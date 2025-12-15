process RIBOCODE_GTFUPDATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocode:1.2.15--pyhfa5458b_0':
        'biocontainers/ribocode:1.2.15--pyhfa5458b_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gtf")                                                  , emit: gtf
    tuple val("${task.process}"), val('ribocode'), eval('RiboCode --version  2>&1') , emit: versions_ribocode, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    GTFupdate \\
        ${gtf} \\
        $args \\
        > ${prefix}_updated.gtf
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_updated.gtf
    """
}
