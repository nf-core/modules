process VIBER {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-cnaqc_r-viber:014077a3164189d5':
        'community.wave.seqera.io/library/r-cnaqc_r-viber:2314592f7d2f9abe'}"

    input:
    tuple val(meta), path(rds_join), val(tumour_samples)

    output:
    tuple val(meta), path("*_viber_best_fit.rds")                , emit: viber_rds
    tuple val(meta), path("*_viber_best_heuristic_fit.rds")      , emit: viber_heuristic_rds
    tuple val(meta), path("*_viber_best_fit_plots.rds")          , emit: viber_plots_rds
    tuple val(meta), path("*_viber_best_heuristic_fit_plots.rds"), emit: viber_heuristic_plots_rds
    tuple val(meta), path("*_viber_report.rds")                  , emit: viber_report_rds
    tuple val(meta), path("*_viber_report.pdf")                  , emit: viber_report_pdf
    tuple val(meta), path("*_viber_report.png")                  , emit: viber_report_png
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "viber_main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_viber_best_fit.rds
    touch ${prefix}_viber_best_heuristic_fit.rds
    touch ${prefix}_viber_best_fit_plots.rds
    touch ${prefix}_viber_best_heuristic_fit_plots.rds
    touch ${prefix}_viber_report.rds
    touch ${prefix}_viber_report.pdf
    touch ${prefix}_viber_report.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VIBER: \$(Rscript -e 'library(VIBER); sessionInfo()\$otherPkgs\$VIBER\$Version')
        cli: \$(Rscript -e 'library(cli); sessionInfo()\$otherPkgs\$cli\$Version')
        dplyr: \$(Rscript -e 'library(dplyr); sessionInfo()\$otherPkgs\$dplyr\$Version')
        ggplot2: \$(Rscript -e 'library(ggplot2); sessionInfo()\$otherPkgs\$ggplot2\$Version')
    END_VERSIONS
    """
}
