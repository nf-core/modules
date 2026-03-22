process UNZIP {
    tag "$archive"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}/"), emit: unzipped_archive
    tuple val("${task.process}"), val('7za'), eval('7za --help | sed -n "s/.*p7zip Version //p" | sed "s/ .*//"'), topic: versions, emit: versions_7za

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if ( archive instanceof List && archive.name.size > 1 ) { error "[UNZIP] error: 7za only accepts a single archive as input. Please check module input." }
    prefix = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.baseName)
    """
    7za \\
        x \\
        -o"${prefix}"/ \\
        $args \\
        $archive
    """

    stub:
    if ( archive instanceof List && archive.name.size > 1 ) { error "[UNZIP] error: 7za only accepts a single archive as input. Please check module input." }
    prefix = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.baseName)
    """
    mkdir "${prefix}"
    """
}
