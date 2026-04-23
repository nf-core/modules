process KMCP_PROFILE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmcp:0.9.4--h9ee0642_0'
        : 'quay.io/biocontainers/kmcp:0.9.4--h9ee0642_0'}"

    input:
    tuple val(meta), path(search_results)
    path db

    output:
    tuple val(meta), path("*.profile"), emit: profile
    tuple val("${task.process}"), val('kmcp'), eval("kmcp version 2>&1 | sed 's/^.*kmcp v//'"), emit: versions_kmcp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    taxid=`find -L ${db} -name "*map"`
    taxdump=`find -L ${db}/*/ -type d  -not -name "R001"`
    kmcp \\
        profile \\
        ${args} \\
        -X \$taxdump \\
        -T \$taxid \\
        -j ${task.cpus} \\
        -o ${prefix}.profile \\
        ${search_results}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile
    """
}
