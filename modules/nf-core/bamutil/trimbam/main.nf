process BAMUTIL_TRIMBAM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamutil:1.0.15--h2e03b76_1' :
        'biocontainers/bamutil:1.0.15--h2e03b76_1' }"

    input:
    tuple val(meta), path(bam), val(trim_left), val(trim_right)

    output:
    tuple val(meta), path("*.bam")                                                                       , emit: bam
    tuple val("${task.process}"), val('bamutil'), eval("bam trimBam 2>&1 | head -1 | sed 's/^Version: //;s/;.*//'"), emit: versions_bamutil, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam \\
        trimBam \\
        $bam \\
        ${prefix}.bam \\
        $args \\
        -L $trim_left \\
        -R $trim_right
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_trimbam"
    """
    touch ${prefix}.bam
    """
}
