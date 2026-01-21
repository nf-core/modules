process SAMTOOLS_QUICKCHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23--h96c455f_0':
        'biocontainers/samtools:1.23--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), env("EXIT_CODE"),   emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools quickcheck \\
        $args \\
        $bam || EXIT_CODE=\$?

    EXIT_CODE=\${EXIT_CODE:-0}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    EXIT_CODE=0
    echo $args
    """
}
