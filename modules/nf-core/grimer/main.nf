process GRIMER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grimer:1.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/grimer:1.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_table)
    path sample_metadata
    path config

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val("${task.process}"), val('grimer'), eval('grimer --version 2>&1 | tail -1'), emit: versions_grimer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def m_arg  = sample_metadata ? "-m ${sample_metadata}" : ''
    def c_arg  = config          ? "-c ${config}"   : ''
    """
    grimer \\
        -i ${input_table} \\
        ${m_arg} \\
        ${c_arg} \\
        -o ${prefix}.html \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    """
}
