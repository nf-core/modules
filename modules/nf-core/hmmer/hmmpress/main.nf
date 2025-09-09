process HMMER_HMMPRESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'biocontainers/hmmer:3.4--hdbdd923_1' }"

    input:
    tuple val(meta), path(hmmfile)

    output:
    tuple val(meta), path("*.h3?"), emit: compressed_db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    hmmpress \\
        $args \\
        ${hmmfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(echo \$(hmmpress -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "stub"

    """
    touch ${prefix}.h3m
    touch ${prefix}.h3i
    touch ${prefix}.h3f
    touch ${prefix}.h3p

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(echo \$(hmmpress -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
    END_VERSIONS
    """
}
