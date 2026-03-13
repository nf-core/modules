process OPENMS_TEXTEXPORTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0'
        : 'biocontainers/openms:3.5.0--h78fb946_0'}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    TextExporter \\
        -in ${input_file} \\
        -out ${prefix}.tsv \\
        -threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
