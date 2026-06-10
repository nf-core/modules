process DEEPMASED_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deepmased:0.3.1--pyh5ca1d4c_0':
        'quay.io/biocontainers/deepmased:0.3.1--pyh5ca1d4c_0' }"

    input:
    tuple val(meta), path(feature_file_table), path(feature_files)

    output:
    tuple val(meta), path("*_predictions.tsv"), emit: predictions
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    DeepMAsED predict \\
        ${feature_file_table} \\
        --n-procs ${task.cpus} \\
        --save-name ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepmased: $VERSION
        setuptools: 78.1
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.1'
    """
    touch ${prefix}_predictions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepmased: $VERSION
        setuptools: 78.1
    END_VERSIONS
    """
}
