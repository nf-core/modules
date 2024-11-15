process NONPAREIL_NONPAREILCURVESR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.5.5--r43hdcf5f25_0':
        'biocontainers/nonpareil:3.5.5--r43hdcf5f25_0' }"

    input:
    tuple val(meta), path(npos)

    output:
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.tsv" ), emit: tsv , optional: true
    tuple val(meta), path("*.csv" ), emit: csv , optional: true
    tuple val(meta), path("*.pdf" ), emit: pdf , optional: true

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    NonpareilCurves.R \\
        $args \\
        --json ${prefix}.json \\
        --tsv ${prefix}.tsv \\
        --csv ${prefix}.csv \\
        --pdf ${prefix}.pdf \\
        $npos

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(echo \$(nonpareil -V 2>&1) | sed 's/Nonpareil v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    touch ${prefix}.tsv
    touch ${prefix}.csv
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(echo \$(nonpareil -V 2>&1) | sed 's/Nonpareil v//' )
    END_VERSIONS
    """
}
