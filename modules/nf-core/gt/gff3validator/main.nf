process GT_GFF3VALIDATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path('*.success.log')  , emit: success_log , optional: true
    tuple val(meta), path('*.error.log')    , emit: error_log   , optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gt \\
        gff3validator \\
        "$gff3" \\
        > "${prefix}.stdout" \\
        2> >(tee "${prefix}.stderr" >&2) \\
        || echo "Errors from gt-gff3validator printed to ${prefix}.error.log"

    if grep -q "input is valid GFF3" "${prefix}.stdout"; then
        echo "Validation successful..."
        # emit stdout to the success output channel
        mv \\
            "${prefix}.stdout" \\
            "${prefix}.success.log"
    else
        echo "Validation failed..."
        # emit stderr to the error output channel
        mv \\
            "${prefix}.stderr" \\
            "${prefix}.error.log"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.success.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """
}
