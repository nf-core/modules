
process SYLPHTAX_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph-tax:1.1.2--pyhdfd78af_0':
        'biocontainers/sylph-tax:1.1.2--pyhdfd78af_0' }"


    input:
    tuple val(meta), path(sylphtax_reports)
    val data_type

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sylph-tax \\
        merge \\
        ${sylphtax_reports} \\
        --column $data_type \\
        --output ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version 2>&1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p' | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version 2>&1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p' | head -n 1)
    END_VERSIONS
    """
}
