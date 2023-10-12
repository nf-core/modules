process ZIP {
    tag "$prefix"
    label 'process_single'

    conda "conda-forge::p7zip=16.02"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(files, stageAs: "inputs/*")

    output:
    tuple val(meta), path("${prefix}.zip"), emit: zipped_archive
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    prefix = task.ext.prefix ?: ( meta.id ? "${meta.id}" : 'zipped_files')
    """
    7z \\
        a \\
        -l \\
        $args \\
        "${prefix}.zip" ./inputs/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
