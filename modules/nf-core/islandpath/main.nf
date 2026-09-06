process ISLANDPATH {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/islandpath:1.0.6--hdfd78af_0':
        'quay.io/biocontainers/islandpath:1.0.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("*.gff")        , emit: gff
    path "Dimob.log"                      , emit: log
    tuple val("${task.process}"), val('islandpath'), val('1.0.6'), topic: versions, emit: versions_islandpath

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    islandpath \\
        $genome \\
        ${prefix}.gff \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff
    touch Dimob.log
    """
}
