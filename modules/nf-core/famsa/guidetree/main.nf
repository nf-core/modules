
process FAMSA_GUIDETREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/famsa:2.4.1--h9ee0642_0':
        'quay.io/biocontainers/famsa:2.4.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.dnd"), emit: tree
    tuple val("${task.process}"), val('famsa'), eval("famsa -help 2>&1 | sed '2!d;s/.*version //;s/ .*//'"), topic: versions, emit: versions_famsa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    famsa -gt_export \\
        $args \\
        -t ${task.cpus} \\
        ${fasta} \\
        ${prefix}.dnd
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dnd
    """
}
