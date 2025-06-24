process CTREE {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ba99151f85012475e44f3f1d2537912d5b69cf40:2b2850df76b35c8b35e75cf13dce31428d2977d4-0':
        'biocontainers/mulled-v2-ba99151f85012475e44f3f1d2537912d5b69cf40:2b2850df76b35c8b35e75cf13dce31428d2977d4-0' }"

    input:
        tuple val(meta), path(ctree_input)

    output:
        tuple val(meta), path("**ctree_{mobster,VIBER,pyclonevi}.rds")             , emit: ctree_rds       , optional: true
        tuple val(meta), path("**ctree_{mobster,VIBER,pyclonevi}_plots.rds")       , emit: ctree_plots_rds , optional: true
        tuple val(meta), path("**REPORT_plots_ctree_{mobster,VIBER,pyclonevi}.rds"), emit: ctree_report_rds, optional: true
        tuple val(meta), path("**REPORT_plots_ctree_{mobster,VIBER,pyclonevi}.pdf"), emit: ctree_report_pdf, optional: true
        tuple val(meta), path("**REPORT_plots_ctree_{mobster,VIBER,pyclonevi}.png"), emit: ctree_report_png, optional: true
        path "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}_ctree_pyclonevi.rds
    touch ${prefix}_ctree_pyclonevi_plots.rds
    touch ${prefix}_REPORT_plots_ctree_pyclonevi.rds
    touch ${prefix}_REPORT_plots_ctree_pyclonevi.pdf
    touch ${prefix}_REPORT_plots_ctree_pyclonevi.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctree: \$(Rscript -e 'library(ctree); sessionInfo()\$otherPkgs\$ctree\$Version')
        mobster: \$(Rscript -e 'library(mobster); sessionInfo()\$otherPkgs\$mobster\$Version')
        VIBER: \$(Rscript -e 'library(VIBER); sessionInfo()\$otherPkgs\$VIBER\$Version')
        cli: \$(Rscript -e 'library(cli); sessionInfo()\$otherPkgs\$cli\$Version')
        dplyr: \$(Rscript -e 'library(dplyr); sessionInfo()\$otherPkgs\$dplyr\$Version')
        ggplot2: \$(Rscript -e 'library(ggplot2); sessionInfo()\$otherPkgs\$ggplot2\$Version')
    END_VERSIONS
    """
}
