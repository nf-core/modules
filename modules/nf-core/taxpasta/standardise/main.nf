process TAXPASTA_STANDARDISE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::taxpasta=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(profile)
    path taxonomy

    output:
    tuple val(meta), path("*."), emit: profiles
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    taxpasta standardise \\
        $args \\
        '$profile'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxpasta: \$(echo \$(taxpasta --version))
    END_VERSIONS
    """
}
