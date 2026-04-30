process SVDSS_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdss:2.1.1--he17396a_0' :
        'quay.io/biocontainers/svdss:2.1.1--he17396a_0' }"

    input:
    tuple val(meta), path(fasta), path(existing_index)
    val output_format

    output:
    tuple val(meta), path("${prefix}.${output_format}"), emit: index
    tuple val("${task.process}"), val('svdss'), eval("SVDSS --version 2>&1 | sed 's/SVDSS, v//'"), emit: versions_svdss, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def format_flag = output_format == 'fmd' ? '-d' : '-b'
    def existing_index_arg = existing_index ? "-i ${existing_index}" : ''
    """
    SVDSS index \\
        -t ${task.cpus} \\
        ${existing_index_arg} \\
        ${format_flag} ${fasta} \\
        ${args} \\
        > ${prefix}.${output_format}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.${output_format}
    """
}
