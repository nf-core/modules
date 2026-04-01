process METHURATOR_PLOT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methurator:2.1.1--pyhdfd78af_0' :
        'biocontainers/methurator:2.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(summary_report)

    output:
    tuple val(meta), path("plots/*.html")  , emit: plots
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    methurator plot \\
        --summary $summary_report \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methurator: "\$(methurator --version 2>&1 | sed -E 's/.*version[[:space:]]+([0-9.]+).*/\\1/')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "plots/${meta.id}.html"
    """
    mkdir plots/
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methurator: "\$(methurator --version 2>&1 | sed -E 's/.*version[[:space:]]+([0-9.]+).*/\\1/')"
    END_VERSIONS
    """
}
