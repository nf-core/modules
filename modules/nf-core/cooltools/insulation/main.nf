process COOLTOOLS_INSULATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.7.1--py39hff726c5_2' :
        'biocontainers/cooltools:0.7.1--py39hff726c5_2' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*tsv"), emit:tsv
    tuple val(meta), path("*.bw"), emit: bigwig, optional: true
    tuple val("${task.process}"), val('cooltools'), eval("cooltools --version | sed -n 's/cooltools, version //p'"), topic: versions, emit: versions_cooltools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools insulation \\
        -p ${task.cpus} \\
        -o ${prefix}_insulation.tsv \\
        ${cool} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_insulation.tsv
    """

}
