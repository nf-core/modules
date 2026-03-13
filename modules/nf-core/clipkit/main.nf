process CLIPKIT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clipkit:2.11.4--pyhdfd78af_0':
        'biocontainers/clipkit:2.11.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(aln)
    val out_format
    path auxiliary_file

    output:
    tuple val(meta), path("${prefix}.${out_extension}")            , emit: clipkit
    tuple val(meta), path("${prefix}.${out_extension}.log")        , emit: log
    tuple val(meta), path("${prefix}.${out_extension}.txt")        , emit: metadata
    tuple val(meta), path("${prefix}.${out_extension}.complement") , emit: complementary, optional: true
    tuple val(meta), path("${prefix}.${out_extension}.report.json"), emit: json         , optional: true
    tuple val("${task.process}"), val('clipkit'), eval("clipkit --version 2>&1 | sed 's/clipkit //'"), topic: versions, emit: versions_clipkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "clipkit"
    def aux_flag  = auxiliary_file ? "--auxiliary_file ${auxiliary_file}" : ''
    """
    clipkit \\
        $args \\
        $aux_flag \\
        --threads ${task.cpus} \\
        $aln \\
        --output ${prefix}.${out_extension} \\
        --log \\
        2>&1 | tee ${prefix}.${out_extension}.txt
    """

    stub:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "clipkit"
    """
    echo $args

    touch ${prefix}.${out_extension}
    touch ${prefix}.${out_extension}.log
    touch ${prefix}.${out_extension}.txt
    ${args.contains('--complement')  ? "touch ${prefix}.${out_extension}.complement"   : ''}
    ${args.contains('--report_json') ? "touch ${prefix}.${out_extension}.report.json"  : ''}
    """
}
