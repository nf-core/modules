process MOBSTER {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ba99151f85012475e44f3f1d2537912d5b69cf40:2b2850df76b35c8b35e75cf13dce31428d2977d4-0':
        'biocontainers/mulled-v2-ba99151f85012475e44f3f1d2537912d5b69cf40:2b2850df76b35c8b35e75cf13dce31428d2977d4-0' }"

    input:
    tuple val(meta), path(rds_join)

    output:
    tuple val(meta), path("*_mobster_st_fit.rds")      , emit: mobster_rds
    tuple val(meta), path("*_mobster_st_best_fit.rds") , emit: mobster_best_rds
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
    touch ${prefix}_REPORT_plots_mobster.ng

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNAqc: \$(Rscript -e 'library(CNAqc); sessionInfo()\$otherPkgs\$CNAqc\$Version')
        mobster: \$(Rscript -e 'library(mobster); sessionInfo()\$otherPkgs\$mobster\$Version')
    END_VERSIONS
    """
}
