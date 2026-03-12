process SAMTOOLS_CRAMSIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), path("*.size"), emit: size
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    """
    samtools \\
        cram-size \\
        $args \\
        -o ${prefix}.size \\
        $cram
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"
    """
    touch ${prefix}.size
    """
}
