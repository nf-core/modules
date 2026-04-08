process FALINT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fa-lint:1.2.0--he881be0_0':
        'biocontainers/fa-lint:1.2.0--he881be0_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.success.log')  , emit: success_log , optional: true
    tuple val(meta), path('*.error.log')    , emit: error_log   , optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    fa-lint \\
        -threads ${task.cpus} \\
        $args \\
        -fasta $fasta \\
        > >(tee ${prefix}.success.log >&1) \\
        2> >(tee ${prefix}.error.log >&2) \\
        || echo "Errors from fa-lint printed to ${prefix}.error.log"

    if [ \$(cat ${prefix}.error.log | wc -l) -gt 0 ]; then
        echo "Validation failed..."

        rm ${prefix}.success.log
    else
        echo "Validation successful..."

        rm ${prefix}.error.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fa-lint: \$(fa-lint -version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Fasta is valid: $fasta" \\
        > "${prefix}.success.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fa-lint: \$(fa-lint -version)
    END_VERSIONS
    """
}
