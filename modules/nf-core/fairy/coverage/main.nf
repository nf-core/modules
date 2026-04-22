process FAIRY_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fairy:0.5.8--hc1c3326_0':
        'quay.io/biocontainers/fairy:0.5.8--hc1c3326_0' }"

    input:
    tuple val(meta), path(sketches), path(contigs)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: coverage
    tuple val("${task.process}"), val('fairy'), eval("fairy --version 2>&1 | sed 's/fairy //'"), topic: versions, emit: versions_fairy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fairy coverage \\
        ${sketches} \\
        ${contigs} \\
        -t ${task.cpus} \\
        -o ${prefix}.tsv \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "$args"

    touch ${prefix}.tsv
    """
}
