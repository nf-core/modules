process ABRICATE_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate%3A1.0.1--ha8f3691_1':
        'biocontainers/abricate:1.0.1--ha8f3691_1' }"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val("${task.process}"), val('abricate'), eval("echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' "), emit: versions_abricate, topic: versions

    tuple val(meta), path("*.txt"), emit: report
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    abricate \\
        --summary \\
        ${reports} > ${prefix}.txt


    """
}
