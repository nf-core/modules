process MOBSTER {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-cnaqc_r-mobster_r-cli_r-dplyr_r-ggplot2:9f1d68529fc936de':
        'community.wave.seqera.io/library/r-cnaqc_r-mobster_r-cli_r-dplyr_r-ggplot2:96c0dbada588b39a' }"

    input:
    tuple val(meta), path(rds_join)

    output:
    tuple val(meta), path("*_mobster_st_fit.rds")       , emit: mobster_rds
    tuple val(meta), path("*_mobster_st_best_fit.rds")  , emit: mobster_best_rds
    tuple val(meta), path("*_plots.rds")                , emit: mobster_plots_rds
    tuple val(meta), path("*_REPORT_plots_mobster.rds") , emit: mobster_report_rds
    tuple val(meta), path("*_REPORT_plots_mobster.pdf") , emit: mobster_report_pdf
    tuple val(meta), path("*_REPORT_plots_mobster.png") , emit: mobster_report_png
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}_mobster_st_fit.rds
    touch ${prefix}_mobster_st_best_fit.rds
    touch ${prefix}_plots.rds
    touch ${prefix}_REPORT_plots_mobster.rds
    touch ${prefix}_REPORT_plots_mobster.pdf
    touch ${prefix}_REPORT_plots_mobster.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNAqc: \$(Rscript -e 'library(CNAqc); sessionInfo()\$otherPkgs\$CNAqc\$Version')
        mobster: \$(Rscript -e 'library(mobster); sessionInfo()\$otherPkgs\$mobster\$Version')
        cli: \$(Rscript -e 'library(cli); sessionInfo()\$otherPkgs\$cli\$Version')
        dplyr: \$(Rscript -e 'library(dplyr); sessionInfo()\$otherPkgs\$dplyr\$Version')
        ggplot2: \$(Rscript -e 'library(ggplot2); sessionInfo()\$otherPkgs\$ggplot2\$Version')
    END_VERSIONS
    """
}
