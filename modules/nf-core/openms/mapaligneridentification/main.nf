process OPENMS_MAPALIGNERIDENTIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(id_files)

    output:
    tuple val(meta), path("*.trafoXML"), emit: trafoxml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def trafo_out = id_files.collect { "${it.baseName}.trafoXML" }.join(' ')
    """
    MapAlignerIdentification \\
        -in $id_files \\
        -trafo_out $trafo_out \\
        -threads $task.cpus \\
        $args
    """

    stub:
    def trafo_out = id_files.collect { "${it.baseName}.trafoXML" }.join(' ')
    """
    touch $trafo_out
    """
}
