process FASTAVALIDATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/py_fasta_validator:0.6--py37h595c7a6_0':
        'biocontainers/py_fasta_validator:0.6--py37h595c7a6_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.success.log")  , emit: success_log , optional: true
    tuple val(meta), path("*.error.log")    , emit: error_log   , optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    py_fasta_validator \\
        -f $fasta \\
        2> "${prefix}.error.log" \\
        || echo "Errors from fasta_validate printed to ${prefix}.error.log"

    if [ \$(cat "${prefix}.error.log" | wc -l) -gt 0 ]; then
        echo "Validation failed..."

        cat \\
            "${prefix}.error.log"
    else
        echo "Validation successful..."

        mv \\
            "${prefix}.error.log" \\
            fasta_validate.stderr

        echo "Validation successful..." \\
            > "${prefix}.success.log"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        py_fasta_validator: \$(py_fasta_validator -v | sed 's/.* version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Validation successful..." \\
        > "${prefix}.success.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        py_fasta_validator: \$(py_fasta_validator -v | sed 's/.* version //')
    END_VERSIONS
    """
}
