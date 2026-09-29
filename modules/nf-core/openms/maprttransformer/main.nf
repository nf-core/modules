process OPENMS_MAPRTTRANSFORMER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(in_file), path(trafoxml)

    output:
    tuple val(meta), path("${prefix}.${in_file.extension}"), emit: aligned
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_aligned"
    if ("$in_file" == "${prefix}.${in_file.extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    MapRTTransformer \\
        -in $in_file \\
        -trafo_in $trafoxml \\
        -out ${prefix}.${in_file.extension} \\
        -threads $task.cpus \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_aligned"
    if ("$in_file" == "${prefix}.${in_file.extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch ${prefix}.${in_file.extension}
    """
}
