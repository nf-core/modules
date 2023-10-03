process HAMRONIZATION_RGI {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::hamronization=1.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hamronization:1.1.1--pyhdfd78af_0':
        'biocontainers/hamronization:1.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(report)
    val(format)
    val(software_version)
    val(reference_db_version)

    output:
    tuple val(meta), path("*.json") , optional: true, emit: json
    tuple val(meta), path("*.tsv")  , optional: true, emit: tsv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hamronize \\
        rgi \\
        ${report} \\
        $args \\
        --format ${format} \\
        --analysis_software_version ${software_version} \\
        --reference_database_version ${reference_db_version} \\
        --input_file_name ${prefix} \\
        > ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: \$(echo \$(hamronize --version 2>&1) | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
