process FAIRY_SKETCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fairy:0.5.8--hc1c3326_0' :
        'quay.io/biocontainers/fairy:0.5.8--hc1c3326_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}/*.bcsp"), emit: sketch
    tuple val("${task.process}"), val("fairy"), eval("fairy --version 2>&1 | sed 's/fairy //'"), topic: versions, emit: versions_fairy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ?
        "-r ${reads}" :
        "-1 ${reads[0]} -2 ${reads[1]}"
    """
    mkdir -p ${prefix}

    fairy \\
        sketch \\
        ${input_reads} \\
        -d ${prefix} \\
        $args
    """

    stub:
    prefix     = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    """
    echo "$args"
    mkdir -p ${prefix}
    touch ${prefix}/${meta.id}.bcsp
    """
}
