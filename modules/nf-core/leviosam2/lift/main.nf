process LEVIOSAM2_LIFT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leviosam2:0.4.2--h4ac6f70_0':
        'quay.io/biocontainers/leviosam2:0.4.2--h4ac6f70_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(clft)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('leviosam2'), eval("leviosam2 --version"), emit: versions_leviosam2, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    leviosam2 \\
        lift \\
        -t ${task.cpus} \\
        -C ${clft} \\
        -a ${input} \\
        -p ${prefix} \\
        -O bam \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
