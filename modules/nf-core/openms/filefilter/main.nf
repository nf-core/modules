process OPENMS_FILEFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.mzML"),         emit: mzml, optional: true
    tuple val(meta), path("*.featureXML"),   emit: featurexml, optional: true
    tuple val(meta), path("*.consensusXML"), emit: consensusxml, optional: true
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "${file.getExtension()}"
    if ("$file" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    FileFilter \\
        -in $file \\
        -out ${prefix}.${suffix} \\
        -threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "${file.getExtension()}"
    if ("$file" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    touch ${prefix}.${suffix}
    """
}
