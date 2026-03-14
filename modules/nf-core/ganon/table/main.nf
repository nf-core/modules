process GANON_TABLE {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ganon:2.1.0--py310hab1bfa5_1'
        : 'biocontainers/ganon:2.1.0--py310hab1bfa5_1'}"

    input:
    tuple val(meta), path(tre)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('ganon'), eval("ganon --version 2>1 | sed 's/.*ganon //g'"), emit: versions_ganon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ganon \\
        table \\
        --input ${tre} \\
        --output-file ${prefix}.txt \\
        ${args}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    """
}
