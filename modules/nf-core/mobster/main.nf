process MOBSTER {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9ef4992af6754a36358a4fc7a52bdc3928c012070f06140a44ddcf174da4e62/data':
        'community.wave.seqera.io/library/r-cnaqc_r-mobster_r-cli_r-dplyr_r-ggplot2:5a94f700b38065ea' }"

    input:
    tuple val(meta), path(rds_join)

    output:
    tuple val(meta), path("*_mobster_fit.rds")        , emit: mobster_rds
    tuple val(meta), path("*_mobster_best_fit.rds")   , emit: mobster_best_rds
    tuple val(meta), path("*_mobster_best_plots.rds") , emit: mobster_best_plots_rds
    tuple val(meta), path("*_mobster_report.rds")     , emit: mobster_report_rds
    tuple val(meta), path("*_mobster_report.pdf")     , emit: mobster_report_pdf
    tuple val(meta), path("*_mobster_report.png")     , emit: mobster_report_png
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}_mobster_fit.rds
    touch ${prefix}_mobster_best_fit.rds
    touch ${prefix}_mobster_best_plots.rds
    touch ${prefix}_mobster_report.rds
    touch ${prefix}_mobster_report.pdf
    touch ${prefix}_mobster_report.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNAqc: \$(Rscript -e "library(CNAqc); cat(as.character(packageVersion('CNAqc')))")
        mobster: \$(Rscript -e "library(mobster); cat(as.character(packageVersion('mobster')))")
        cli: \$(Rscript -e "library(cli); cat(as.character(packageVersion('cli')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}
