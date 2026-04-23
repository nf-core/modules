process KMCP_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmcp:0.9.4--h9ee0642_0'
        : 'quay.io/biocontainers/kmcp:0.9.4--h9ee0642_0'}"

    input:
    tuple val(meta), path(search_out)

    output:
    tuple val(meta), path("*.gz"), emit: result
    tuple val("${task.process}"), val('kmcp'), eval("kmcp version 2>&1 | sed 's/^.*kmcp v//'"), emit: versions_kmcp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.merged"
    """
    kmcp \\
        merge \\
        ${args} \\
        --threads ${task.cpus} \\
        --out-file ${prefix}.gz \\
        ${search_out}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.merged"
    """
    echo "" | gzip > ${prefix}.gz
    """
}
