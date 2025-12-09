
process CLUSTY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clusty:1.2.2--h9ee0642_0':
        'biocontainers/clusty:1.2.2--h9ee0642_0' }"

    input:
    tuple val(meta), path(distances)
    tuple val(meta2), path(objects)

    output:
    tuple val(meta), path("*.tsv"), emit: assignments
    tuple val("${task.process}"), val('clusty'), eval('echo $(clusty --version 2>&1)'), topic: versions, emit: versions_clusty

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def objects_arg = objects ? "--objects-file $objects" : ""

    if ("${distances}" == "${prefix}.tsv") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    clusty \\
        $args \\
        -t $task.cpus \\
        ${objects_arg} \\
        ${distances} \\
        ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def objects_arg = objects ? "--objects-file $objects" : ""

    if ("${distances}" == "${prefix}.tsv") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    echo $args
    echo ${objects_arg}
    touch ${prefix}.tsv
    """
}
