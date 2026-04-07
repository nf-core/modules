process CTREE {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/4084291fbed2d8371cc9dd8c53a422d0731a27b2366d0dc0069c0fc0fac314bb/data':
        'community.wave.seqera.io/library/r-ctree_r-mobster_r-viber_r-cli_pruned:48299db4104e296b' }"

    input:
    tuple val(meta), path(ctree_input)

    output:
    tuple val(meta), path("**ctree_{mobster,VIBER,pyclonevi}.rds")       , emit: ctree_rds       , optional: true
    tuple val(meta), path("**ctree_{mobster,VIBER,pyclonevi}_plots.rds") , emit: ctree_plots_rds , optional: true
    tuple val(meta), path("**ctree_{mobster,VIBER,pyclonevi}_report.rds"), emit: ctree_report_rds, optional: true
    tuple val(meta), path("**ctree_{mobster,VIBER,pyclonevi}_report.pdf"), emit: ctree_report_pdf, optional: true
    tuple val(meta), path("**ctree_{mobster,VIBER,pyclonevi}_report.png"), emit: ctree_report_png, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    # outputs for mobster rds input
    mkdir stub/
    touch stub/${prefix}_ctree_mobster.rds
    touch stub/${prefix}_ctree_mobster_plots.rds
    touch stub/${prefix}_ctree_mobster_report.rds
    touch stub/${prefix}_ctree_mobster_report.pdf
    touch stub/${prefix}_ctree_mobster_report.png

    # outputs for pyclonevi tsv input
    touch ${prefix}_ctree_pyclonevi.rds
    touch ${prefix}_ctree_pyclonevi_plots.rds
    touch ${prefix}_ctree_pyclonevi_report.rds
    touch ${prefix}_ctree_pyclonevi_report.pdf
    touch ${prefix}_ctree_pyclonevi_report.png

    # outputs for VIBER rds input
    touch ${prefix}_ctree_VIBER.rds
    touch ${prefix}_ctree_VIBER_plots.rds
    touch ${prefix}_ctree_VIBER_report.pdf
    touch ${prefix}_ctree_VIBER_report.png
    touch ${prefix}_ctree_VIBER_report.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctree: \$(Rscript -e "library(ctree); cat(as.character(packageVersion('ctree')))")
        mobster: \$(Rscript -e "library(mobster); cat(as.character(packageVersion('mobster')))")
        VIBER: \$(Rscript -e "library(VIBER); cat(as.character(packageVersion('VIBER')))")
        cli: \$(Rscript -e "library(cli); cat(as.character(packageVersion('cli')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        ggpubr: \$(Rscript -e "library(ggpubr); cat(as.character(packageVersion('ggpubr')))")
    END_VERSIONS
    """
}
