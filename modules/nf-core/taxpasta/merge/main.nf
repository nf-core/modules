process TAXPASTA_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxpasta:0.7.0--pyhdfd78af_0':
        'biocontainers/taxpasta:0.7.0--pyhdfd78af_0' }"


    input:
    tuple val(meta), path(profiles)
    val profiler
    val format
    path taxonomy
    path samplesheet

    output:
    tuple val(meta), path("*.{tsv,csv,arrow,parquet,biom}"), emit: merged_profiles
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // N.B.: Taxpasta requires a --output
    // option with the output file name. The desired format will be parsed from
    // the name and should correspond to the output pattern specified above,
    // e.g., "--output ${task.ext.prefix}.tsv".
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxonomy_option = taxonomy ? "--taxonomy ${taxonomy}" : ''
    def samplesheet_input = samplesheet ? "-s ${samplesheet}" : ''
    """
    taxpasta merge \\
        --profiler $profiler \\
        --output ${prefix}.${format} \\
        $args \\
        $taxonomy_option \\
        $samplesheet_input \\
        $profiles


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxpasta: \$(taxpasta --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxonomy_option = taxonomy ? "--taxonomy ${taxonomy}" : ''
    def samplesheet_input = samplesheet ? "-s ${samplesheet}" : ''
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxpasta: \$(taxpasta --version)
    END_VERSIONS
    """
}
