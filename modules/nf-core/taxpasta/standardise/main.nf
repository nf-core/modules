process TAXPASTA_STANDARDISE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::taxpasta=0.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxpasta:0.6.0--pyhdfd78af_0':
        'biocontainers/taxpasta:0.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(profile)
    path taxonomy

    output:
    tuple val(meta), path("*.{tsv,csv,arrow,parquet,biom}"), emit: standardised_profile
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // N.B.: Taxpasta requires a --profiler option and will fail without it.
    // This must be specified via a `nextflow.config` or `modules.config`, for
    // example, as "--profiler kraken2". Additionally, it requires a --output
    // option with the output file name. The desired format will be parsed from
    // the name and should correspond to the output pattern specified above,
    // e.g., "--output ${task.ext.prefix}.tsv".
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxonomy_option = taxonomy ? "--taxonomy ${taxonomy}" : ''
    """
    taxpasta standardise \\
        $args \\
        $taxonomy_option \\
        $profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxpasta: \$(taxpasta --version)
    END_VERSIONS
    """
}
