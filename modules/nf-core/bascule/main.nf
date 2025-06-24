process BASCULE {
    tag "$meta.id"
    label 'process_low'

    // TODO add bioconda package and containers
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(count_matrices)  // count_matrices = folder with .csv files

    output:
    tuple val(meta), path("*.rds")      , emit: bascule_rds
    tuple val(meta), path("*_plots.rds"), emit: bascule_plots_rds
    tuple val(meta), path("*_plots.pdf"), emit: bascule_plots_pdf
    tuple val(meta), path("*_plots.png"), emit: bascule_plots_png
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.rds
    touch ${prefix}_plots.rds
    touch ${prefix}_plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bascule: \$(bascule --version)
    END_VERSIONS
    """
}
