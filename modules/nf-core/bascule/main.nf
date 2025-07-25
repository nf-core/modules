process BASCULE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-bascule_compilers_python_r-ggplot2_pruned:51342d731abb95e0':
        'community.wave.seqera.io/library/r-bascule_compilers_python_r-ggplot2_pruned:00e043c5dd0fe33f' }"

    input:
    tuple val(meta), path(counts_matrices)  // counts_matrices = folder with .csv files

    output:
    tuple val(meta), path("*_fit_bascule.rds")  , emit: bascule_rds
    tuple val(meta), path("*_plots_bascule.rds"), emit: bascule_plots_rds
    tuple val(meta), path("*_plots_bascule.pdf"), emit: bascule_plots_pdf
    tuple val(meta), path("*_plots_bascule.png"), emit: bascule_plots_png
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export RETICULATE_PYTHON=\$(type -p python3)
    """
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}_fit_bascule.rds
    touch ${prefix}_plots_bascule.rds
    touch ${prefix}_plots_bascule.pdf
    touch ${prefix}_plots_bascule.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bascule: \$(Rscript -e 'library(bascule); sessionInfo()\$otherPkgs\$bascule\$Version')
        ggplot2: \$(Rscript -e 'library(ggplot2); sessionInfo()\$otherPkgs\$ggplot2\$Version')
        reticulate: \$(Rscript -e 'library(reticulate); sessionInfo()\$otherPkgs\$reticulate\$Version')
        tidyverse: \$(Rscript -e 'library(tidyverse); sessionInfo()\$otherPkgs\$tidyverse\$Version')
    END_VERSIONS
    """
}
