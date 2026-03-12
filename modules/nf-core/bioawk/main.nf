process BIOAWK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h5bf99c6_6' :
        'biocontainers/bioawk:1.0--h5bf99c6_6' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz"), emit: output
    tuple val("${task.process}"), val('bioawk'), val(VERSION), emit: versions_bioawk, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${input}" == "${prefix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate."

    def VERSION = '1.0' // Tool does not report version via CLI

    """
    bioawk \
        $args \
        $input \
        > ${prefix}

    gzip ${prefix}
    """
}
