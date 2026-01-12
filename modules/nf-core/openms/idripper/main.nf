process OPENMS_IDRIPPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0':
        'biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(merged_idxml)

    output:
    tuple val(meta), path("*.idXML"), emit: idxmls
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | grep -E '^Version' | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    IDRipper \\
        -in $merged_idxml \\
        -out . \\
        -threads $task.cpus \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_1.idXML
    touch ${prefix}_2.idXML
    """
}
