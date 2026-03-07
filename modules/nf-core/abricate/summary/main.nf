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
    tuple val(meta), path("*.txt"), emit: report
    tuple val("${task.process}"), val('abricate'), eval("abricate --version | sed 's/^.* //' "), emit: versions_abricate, topic: versions

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
