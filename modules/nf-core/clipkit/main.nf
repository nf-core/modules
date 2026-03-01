process CLIPKIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clipkit:2.11.4--pyhdfd78af_0':
        'biocontainers/clipkit:2.11.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(aln)
    val out_format
    path auxiliary_file

    output:
    tuple val(meta), path("${prefix}.${out_extension}")           , emit: clipkit
    tuple val(meta), path("${prefix}.${out_extension}.log")       , emit: log          , optional: true
    tuple val(meta), path("${prefix}.${out_extension}.complement"), emit: complementary, optional: true
    tuple val("${task.process}"), val('clipkit'), eval("clipkit --version 2>&1 | sed 's/clipkit //'"), topic: versions, emit: versions_clipkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "clipkit"
    def aux_flag = auxiliary_file ? "-a ${auxiliary_file}" : ''
    """
    clipkit \\
        $args \\
        $aln \\
        -o ${prefix}.${out_extension} \\
        -l \\
        ${aux_flag}
    """

    stub:
    prefix        = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "clipkit"
    """
    touch ${prefix}.${out_extension}
    touch ${prefix}.${out_extension}.log
    touch ${prefix}.${out_extension}.complement
    """
}
