process PAIRIX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairix:0.3.7--py36h30a8e3e_3' :
        'quay.io/biocontainers/pairix:0.3.7--py36h30a8e3e_3' }"

    input:
    tuple val(meta), path(pair)

    output:
    tuple val(meta), path(pair), path("*.px2"), emit: index
    tuple val("${task.process}"), val('pairix'), eval('echo "$(pairix --help 2>&1)" | sed "s/^.*Version: //; s/Usage.*$//"'), topic: versions, emit: versions_pairix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    pairix \\
        $args \\
        $pair
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${pair}.px2
    """
}
