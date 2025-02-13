
process SYLPHTAX_TAXPROF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph-tax:1.1.2--pyhdfd78af_0':
        'biocontainers/sylph-tax:1.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sylph_results)
    path taxonomy

    output:
    tuple val(meta), path("*.sylphmpa")         , emit: taxprof_output
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sylph-tax \\
        taxprof \\
        $sylph_results \\
        $args \\
        -o ${prefix}.sylphmpa \\
        -t $taxonomy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version 2>&1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p' | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sylphmpa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version 2>&1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p' | head -n 1)
    END_VERSIONS
    """
}
