process TAXPASTA_STANDARDISE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::taxpasta=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(profile)
    val output_format
    path taxonomy

    output:
    tuple val(meta), path("*.{tsv,csv,arrow,parquet,biom}"), emit: profiles
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Taxpasta requires a --profiler option and will fail without it.
    // That needs to be configured since we can't set a default here.
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    taxpasta standardise \\
        $args \\
        --output '${prefix}.${output_format}'
        '$profile'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxpasta: \$(echo \$(taxpasta --version))
    END_VERSIONS
    """
}
