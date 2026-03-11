process ARIA2 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aria2:1.36.0' :
        'biocontainers/aria2:1.36.0' }"

    input:
    tuple val(meta), val(source_url)

    output:
    tuple val(meta), path("$downloaded_file"), emit: downloaded_file
    tuple val("${task.process}"), val("aria2"), eval("aria2c --version 2>&1 | grep 'aria2 version' | cut -f3 -d ' '"), emit: versions_aria2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    downloaded_file = source_url.split("/")[-1]

    """
    aria2c \\
        --check-certificate=false \\
        ${args} \\
        ${source_url}
    """

    stub:
    downloaded_file = source_url.split("/")[-1]

    """
    touch ${downloaded_file}
    """
}
