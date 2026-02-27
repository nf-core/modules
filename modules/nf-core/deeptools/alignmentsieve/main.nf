process DEEPTOOLS_ALIGNMENTSIEVE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0':
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(input), path(input_index)

    output:
    tuple val(meta), path("*_as.bam") , emit: bam
    tuple val("${task.process}"), val('deeptools'), eval('alignmentSieve --version | sed "s/alignmentSieve //g"') , emit: versions_deeptools, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'") , emit: versions_samtools, topic: versions
    path  "*_log.txt"                 , emit: logs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    alignmentSieve \\
        $args \\
        -b $input \\
        -o ${prefix}_as.bam \\
        --filterMetrics ${prefix}_log.txt \\
        --numberOfProcessors $task.cpus
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_as.bam
    touch ${prefix}_log.txt
    """
}
